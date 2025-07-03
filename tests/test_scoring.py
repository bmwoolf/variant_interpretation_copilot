"""
Unit tests for ACMG scoring module.
"""

import pytest

from vcf_copilot.scoring import ACMGClassifier, ACMGCriteriaResult
from vcf_copilot.models import Variant, VariantType, ACMGClassification, ACMGCriteria


class TestACMGClassifier:
    """Test ACMG classifier functionality."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.classifier = ACMGClassifier()
    
    def test_classify_pathogenic_variant(self):
        """Test classification of a pathogenic variant."""
        # Create a pathogenic variant (BRCA1 LOF)
        variant = Variant(
            chrom="chr17",
            pos=43057093,
            ref="A",
            alt="T",
            variant_type=VariantType.SNV,
            gene="BRCA1",
            impact="HIGH",
            hgvs_p="p.Arg1751Ter",
            gnomad_af=0.000001,  # Very rare
            cadd_score=25.0,
            polyphen_score=0.95,
            sift_score=0.01
        )
        
        result = self.classifier.classify(variant)
        
        assert isinstance(result.classification, ACMGClassification)
        assert result.score > 0  # Should have positive score
        assert len(result.criteria) > 0  # Should have applied criteria
        assert result.confidence > 0
    
    def test_classify_benign_variant(self):
        """Test classification of a benign variant."""
        # Create a benign variant (common in population)
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            variant_type=VariantType.SNV,
            gene="GENE1",
            gnomad_af=0.1,  # Common variant
            cadd_score=5.0,
            polyphen_score=0.1,
            sift_score=0.5
        )
        
        result = self.classifier.classify(variant)
        
        assert isinstance(result.classification, ACMGClassification)
        assert result.score < 0  # Should have negative score
        assert result.confidence > 0
    
    def test_evaluate_PVS1(self):
        """Test PVS1 criteria evaluation."""
        # LOF variant in BRCA1 (known LOF mechanism)
        variant = Variant(
            chrom="chr17",
            pos=43057093,
            ref="A",
            alt="T",
            variant_type=VariantType.SNV,
            gene="BRCA1",
            impact="HIGH",
            hgvs_p="p.Arg1751Ter"
        )
        
        result = self.classifier._evaluate_PVS1(variant)
        
        assert isinstance(result, ACMGCriteriaResult)
        assert result.criteria == ACMGCriteria.PVS1
        assert result.met  # Should be met for BRCA1 LOF
        assert result.strength == "Very Strong"
    
    def test_evaluate_PM2(self):
        """Test PM2 criteria evaluation."""
        # Rare variant
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            variant_type=VariantType.SNV,
            gnomad_af=0.00001  # Very rare
        )
        
        result = self.classifier._evaluate_PM2(variant)
        
        assert isinstance(result, ACMGCriteriaResult)
        assert result.criteria == ACMGCriteria.PM2
        assert result.met  # Should be met for rare variant
        assert result.strength == "Moderate"
    
    def test_evaluate_PP3(self):
        """Test PP3 criteria evaluation."""
        # Variant with multiple deleterious predictions
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            variant_type=VariantType.SNV,
            cadd_score=25.0,
            polyphen_score=0.95,
            sift_score=0.01
        )
        
        result = self.classifier._evaluate_PP3(variant)
        
        assert isinstance(result, ACMGCriteriaResult)
        assert result.criteria == ACMGCriteria.PP3
        assert result.met  # Should be met with multiple deleterious predictions
        assert result.strength == "Supporting"
    
    def test_evaluate_BA1(self):
        """Test BA1 criteria evaluation."""
        # Common variant
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            variant_type=VariantType.SNV,
            gnomad_af=0.1  # Common variant
        )
        
        result = self.classifier._evaluate_BA1(variant)
        
        assert isinstance(result, ACMGCriteriaResult)
        assert result.criteria == ACMGCriteria.BA1
        assert result.met  # Should be met for common variant
        assert result.strength == "Standalone"
    
    def test_calculate_score(self):
        """Test ACMG score calculation."""
        # Create mock criteria results
        results = [
            ACMGCriteriaResult(
                criteria=ACMGCriteria.PVS1,
                met=True,
                strength="Very Strong",
                reasoning="Test",
                evidence={}
            ),
            ACMGCriteriaResult(
                criteria=ACMGCriteria.PM2,
                met=True,
                strength="Moderate",
                reasoning="Test",
                evidence={}
            )
        ]
        
        score = self.classifier._calculate_score(results)
        
        assert score == 10  # 8 (PVS1) + 2 (PM2)
    
    def test_determine_classification(self):
        """Test classification determination from score."""
        # Pathogenic
        classification = self.classifier._determine_classification(8)
        assert classification == ACMGClassification.PATHOGENIC
        
        # Likely Pathogenic
        classification = self.classifier._determine_classification(6)
        assert classification == ACMGClassification.LIKELY_PATHOGENIC
        
        # Uncertain Significance
        classification = self.classifier._determine_classification(0)
        assert classification == ACMGClassification.UNCERTAIN_SIGNIFICANCE
        
        # Likely Benign
        classification = self.classifier._determine_classification(-6)
        assert classification == ACMGClassification.LIKELY_BENIGN
        
        # Benign
        classification = self.classifier._determine_classification(-8)
        assert classification == ACMGClassification.BENIGN
    
    def test_is_lof_variant(self):
        """Test LOF variant detection."""
        # High impact variant
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="T",
            variant_type=VariantType.SNV,
            impact="HIGH"
        )
        
        assert self.classifier._is_lof_variant(variant)
        
        # Stop gain in protein
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="T",
            variant_type=VariantType.SNV,
            hgvs_p="p.Arg123Ter"
        )
        
        assert self.classifier._is_lof_variant(variant)
    
    def test_has_lof_mechanism(self):
        """Test LOF mechanism detection."""
        # BRCA1 (known LOF gene)
        variant = Variant(
            chrom="chr17",
            pos=43057093,
            ref="A",
            alt="T",
            variant_type=VariantType.SNV,
            gene="BRCA1"
        )
        
        assert self.classifier._has_lof_mechanism(variant)
        
        # Unknown gene
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="T",
            variant_type=VariantType.SNV,
            gene="UNKNOWN"
        )
        
        assert not self.classifier._has_lof_mechanism(variant) 