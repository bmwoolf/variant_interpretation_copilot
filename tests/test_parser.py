"""
Unit tests for VCF parser module.
"""

import pytest
from pathlib import Path

from vcf_copilot.parser import VCFParser, ValidationResult
from vcf_copilot.models import Variant, VariantType


class TestVCFParser:
    """Test VCF parser functionality."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.parser = VCFParser()
        self.test_vcf_path = Path("data/test.vcf")
    
    def test_parse_vcf_file(self):
        """Test parsing a VCF file."""
        if not self.test_vcf_path.exists():
            pytest.skip("Test VCF file not found")
        
        variants = self.parser.parse(self.test_vcf_path)
        
        assert len(variants) == 3
        assert all(isinstance(v, Variant) for v in variants)
        
        # Check first variant (TP53)
        tp53_variant = variants[0]
        assert tp53_variant.chrom == "chr17"
        assert tp53_variant.pos == 7577120
        assert tp53_variant.ref == "A"
        assert tp53_variant.alt == "G"
        assert tp53_variant.variant_type == VariantType.SNV
        assert tp53_variant.gene == "TP53"
        assert tp53_variant.impact == "HIGH"
        assert tp53_variant.cadd_score == 25.3
    
    def test_validate_vcf_file(self):
        """Test VCF file validation."""
        if not self.test_vcf_path.exists():
            pytest.skip("Test VCF file not found")
        
        result = self.parser.validate(self.test_vcf_path)
        
        assert isinstance(result, ValidationResult)
        assert result.is_valid
        assert result.variant_count == 3
        assert result.sample_count == 1
        assert len(result.errors) == 0
    
    def test_determine_variant_type(self):
        """Test variant type determination."""
        # Test SNV
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="T",
            variant_type=VariantType.SNV
        )
        assert variant.variant_type == VariantType.SNV
        
        # Test INDEL
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="AT",
            variant_type=VariantType.INDEL
        )
        assert variant.variant_type == VariantType.INDEL
    
    def test_parse_info_field(self):
        """Test INFO field parsing."""
        # This would test the _parse_info method
        # Implementation depends on cyvcf2 record structure
        pass
    
    def test_get_header_info(self):
        """Test header information extraction."""
        if not self.test_vcf_path.exists():
            pytest.skip("Test VCF file not found")
        
        header_info = self.parser.get_header_info(self.test_vcf_path)
        
        assert isinstance(header_info, dict)
        assert 'samples' in header_info
        assert 'sample_count' in header_info
        assert header_info['sample_count'] == 1
        assert 'SAMPLE1' in header_info['samples'] 

def test_clinvar_query_chromosome_normalization():
    from vcf_copilot.annotators.clinvar import ClinVarAnnotator
    from vcf_copilot.models import Variant
    annotator = ClinVarAnnotator()
    variant = Variant(chrom='chr17', pos=7577120, ref='A', alt='G', variant_type='SNV')
    query = annotator._build_search_query(variant)
    assert 'chr17' not in query
    assert '17[chr]' in query

def test_gnomad_query_chromosome_normalization():
    from vcf_copilot.annotators.gnomad import GnomADAnnotator
    from vcf_copilot.models import Variant
    annotator = GnomADAnnotator()
    variant = Variant(chrom='chr17', pos=7577120, ref='A', alt='G', variant_type='SNV')
    query = annotator._build_graphql_query(variant)
    assert 'chr17' not in query
    assert '17-7577120-A-G' in query 