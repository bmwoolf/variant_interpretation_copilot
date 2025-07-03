"""
ACMG classification module for clinical variant interpretation.
"""

import logging
from typing import List, Dict, Optional
from dataclasses import dataclass

from .models import Variant, ACMGClassification, ACMGCriteria, ACMGResult

logger = logging.getLogger(__name__)


@dataclass
class ACMGCriteriaResult:
    """Result of ACMG criteria evaluation."""
    criteria: ACMGCriteria
    met: bool
    strength: str  # "Very Strong", "Strong", "Moderate", "Supporting"
    reasoning: str
    evidence: Dict


class ACMGClassifier:
    """ACMG criteria classifier for variant interpretation."""
    
    def __init__(self):
        """Initialize ACMG classifier."""
        self.logger = logging.getLogger(__name__)
        
        # Define criteria strengths
        self.criteria_strengths = {
            # Pathogenic Very Strong
            ACMGCriteria.PVS1: "Very Strong",
            # Pathogenic Strong
            ACMGCriteria.PS1: "Strong",
            ACMGCriteria.PS2: "Strong", 
            ACMGCriteria.PS3: "Strong",
            ACMGCriteria.PS4: "Strong",
            # Pathogenic Moderate
            ACMGCriteria.PM1: "Moderate",
            ACMGCriteria.PM2: "Moderate",
            ACMGCriteria.PM3: "Moderate",
            ACMGCriteria.PM4: "Moderate",
            ACMGCriteria.PM5: "Moderate",
            ACMGCriteria.PM6: "Moderate",
            # Pathogenic Supporting
            ACMGCriteria.PP1: "Supporting",
            ACMGCriteria.PP2: "Supporting",
            ACMGCriteria.PP3: "Supporting",
            ACMGCriteria.PP4: "Supporting",
            ACMGCriteria.PP5: "Supporting",
            # Benign Standalone
            ACMGCriteria.BA1: "Standalone",
            # Benign Strong
            ACMGCriteria.BS1: "Strong",
            ACMGCriteria.BS2: "Strong",
            ACMGCriteria.BS3: "Strong",
            ACMGCriteria.BS4: "Strong",
            # Benign Supporting
            ACMGCriteria.BP1: "Supporting",
            ACMGCriteria.BP2: "Supporting",
            ACMGCriteria.BP3: "Supporting",
            ACMGCriteria.BP4: "Supporting",
            ACMGCriteria.BP5: "Supporting",
            ACMGCriteria.BP6: "Supporting",
            ACMGCriteria.BP7: "Supporting",
        }
    
    def classify(self, variant: Variant) -> ACMGResult:
        """
        Classify variant according to ACMG criteria.
        
        Args:
            variant: Variant to classify
            
        Returns:
            ACMG classification result
        """
        try:
            # Evaluate all ACMG criteria
            criteria_results = self._evaluate_criteria(variant)
            
            # Calculate ACMG score
            score = self._calculate_score(criteria_results)
            
            # Determine classification
            classification = self._determine_classification(score)
            
            # Generate reasoning
            reasoning = self._generate_reasoning(criteria_results, classification)
            
            # Calculate confidence
            confidence = self._calculate_confidence(criteria_results)
            
            return ACMGResult(
                variant=variant,
                classification=classification,
                criteria=[cr.criteria for cr in criteria_results if cr.met],
                score=score,
                reasoning=reasoning,
                confidence=confidence,
            )
            
        except Exception as e:
            self.logger.error(f"ACMG classification failed for {variant.chrom}:{variant.pos}: {e}")
            # Return uncertain significance as fallback
            return ACMGResult(
                variant=variant,
                classification=ACMGClassification.UNCERTAIN_SIGNIFICANCE,
                criteria=[],
                score=0,
                reasoning=f"Classification failed: {e}",
                confidence=0.0,
            )
    
    def _evaluate_criteria(self, variant: Variant) -> List[ACMGCriteriaResult]:
        """
        Evaluate all ACMG criteria for a variant.
        
        Args:
            variant: Variant to evaluate
            
        Returns:
            List of criteria evaluation results
        """
        results = []
        
        # Pathogenic criteria
        results.append(self._evaluate_PVS1(variant))
        results.append(self._evaluate_PS1(variant))
        results.append(self._evaluate_PS2(variant))
        results.append(self._evaluate_PS3(variant))
        results.append(self._evaluate_PS4(variant))
        results.append(self._evaluate_PM1(variant))
        results.append(self._evaluate_PM2(variant))
        results.append(self._evaluate_PM3(variant))
        results.append(self._evaluate_PM4(variant))
        results.append(self._evaluate_PM5(variant))
        results.append(self._evaluate_PM6(variant))
        results.append(self._evaluate_PP1(variant))
        results.append(self._evaluate_PP2(variant))
        results.append(self._evaluate_PP3(variant))
        results.append(self._evaluate_PP4(variant))
        results.append(self._evaluate_PP5(variant))
        
        # Benign criteria
        results.append(self._evaluate_BA1(variant))
        results.append(self._evaluate_BS1(variant))
        results.append(self._evaluate_BS2(variant))
        results.append(self._evaluate_BS3(variant))
        results.append(self._evaluate_BS4(variant))
        results.append(self._evaluate_BP1(variant))
        results.append(self._evaluate_BP2(variant))
        results.append(self._evaluate_BP3(variant))
        results.append(self._evaluate_BP4(variant))
        results.append(self._evaluate_BP5(variant))
        results.append(self._evaluate_BP6(variant))
        results.append(self._evaluate_BP7(variant))
        
        return results
    
    def _evaluate_PVS1(self, variant: Variant) -> ACMGCriteriaResult:
        """Evaluate PVS1: Null variant in gene where LOF is a known mechanism of disease."""
        met = False
        reasoning = "PVS1: Not met - insufficient evidence for LOF mechanism"
        
        # Check if variant is predicted to cause loss of function
        if self._is_lof_variant(variant):
            # Check if gene has known LOF mechanism
            if self._has_lof_mechanism(variant):
                met = True
                reasoning = "PVS1: Met - null variant in gene with known LOF mechanism"
        
        return ACMGCriteriaResult(
            criteria=ACMGCriteria.PVS1,
            met=met,
            strength=self.criteria_strengths[ACMGCriteria.PVS1],
            reasoning=reasoning,
            evidence={}
        )
    
    def _evaluate_PM2(self, variant: Variant) -> ACMGCriteriaResult:
        """Evaluate PM2: Absent from controls in population databases."""
        met = False
        reasoning = "PM2: Not met - variant present in population databases"
        
        # Check gnomAD frequency
        if variant.gnomad_af is not None:
            if variant.gnomad_af < 0.0001:  # < 0.01%
                met = True
                reasoning = f"PM2: Met - rare variant (AF: {variant.gnomad_af:.6f})"
            else:
                reasoning = f"PM2: Not met - common variant (AF: {variant.gnomad_af:.6f})"
        else:
            reasoning = "PM2: Not evaluated - no population frequency data"
        
        return ACMGCriteriaResult(
            criteria=ACMGCriteria.PM2,
            met=met,
            strength=self.criteria_strengths[ACMGCriteria.PM2],
            reasoning=reasoning,
            evidence={'gnomad_af': variant.gnomad_af}
        )
    
    def _evaluate_PP3(self, variant: Variant) -> ACMGCriteriaResult:
        """Evaluate PP3: Multiple lines of computational evidence support deleterious effect."""
        met = False
        reasoning = "PP3: Not met - insufficient computational evidence"
        
        # Count deleterious predictions
        deleterious_count = 0
        total_count = 0
        
        if variant.cadd_score is not None:
            total_count += 1
            if variant.cadd_score > 20:  # CADD threshold
                deleterious_count += 1
        
        if variant.polyphen_score is not None:
            total_count += 1
            if variant.polyphen_score > 0.908:  # PolyPhen probably damaging
                deleterious_count += 1
        
        if variant.sift_score is not None:
            total_count += 1
            if variant.sift_score < 0.05:  # SIFT deleterious
                deleterious_count += 1
        
        # Require at least 2 deleterious predictions
        if deleterious_count >= 2 and total_count >= 2:
            met = True
            reasoning = f"PP3: Met - {deleterious_count}/{total_count} tools predict deleterious"
        
        return ACMGCriteriaResult(
            criteria=ACMGCriteria.PP3,
            met=met,
            strength=self.criteria_strengths[ACMGCriteria.PP3],
            reasoning=reasoning,
            evidence={
                'cadd_score': variant.cadd_score,
                'polyphen_score': variant.polyphen_score,
                'sift_score': variant.sift_score,
                'deleterious_count': deleterious_count,
                'total_count': total_count
            }
        )
    
    def _evaluate_BA1(self, variant: Variant) -> ACMGCriteriaResult:
        """Evaluate BA1: Allele frequency > 5% in population databases."""
        met = False
        reasoning = "BA1: Not met - variant not common in population"
        
        if variant.gnomad_af is not None and variant.gnomad_af > 0.05:
            met = True
            reasoning = f"BA1: Met - common variant (AF: {variant.gnomad_af:.6f})"
        
        return ACMGCriteriaResult(
            criteria=ACMGCriteria.BA1,
            met=met,
            strength=self.criteria_strengths[ACMGCriteria.BA1],
            reasoning=reasoning,
            evidence={'gnomad_af': variant.gnomad_af}
        )
    
    # Placeholder methods for other criteria
    def _evaluate_PS1(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PS1, False, "Strong", "PS1: Not evaluated", {})
    
    def _evaluate_PS2(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PS2, False, "Strong", "PS2: Not evaluated", {})
    
    def _evaluate_PS3(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PS3, False, "Strong", "PS3: Not evaluated", {})
    
    def _evaluate_PS4(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PS4, False, "Strong", "PS4: Not evaluated", {})
    
    def _evaluate_PM1(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PM1, False, "Moderate", "PM1: Not evaluated", {})
    
    def _evaluate_PM3(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PM3, False, "Moderate", "PM3: Not evaluated", {})
    
    def _evaluate_PM4(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PM4, False, "Moderate", "PM4: Not evaluated", {})
    
    def _evaluate_PM5(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PM5, False, "Moderate", "PM5: Not evaluated", {})
    
    def _evaluate_PM6(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PM6, False, "Moderate", "PM6: Not evaluated", {})
    
    def _evaluate_PP1(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PP1, False, "Supporting", "PP1: Not evaluated", {})
    
    def _evaluate_PP2(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PP2, False, "Supporting", "PP2: Not evaluated", {})
    
    def _evaluate_PP4(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PP4, False, "Supporting", "PP4: Not evaluated", {})
    
    def _evaluate_PP5(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.PP5, False, "Supporting", "PP5: Not evaluated", {})
    
    def _evaluate_BS1(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BS1, False, "Strong", "BS1: Not evaluated", {})
    
    def _evaluate_BS2(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BS2, False, "Strong", "BS2: Not evaluated", {})
    
    def _evaluate_BS3(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BS3, False, "Strong", "BS3: Not evaluated", {})
    
    def _evaluate_BS4(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BS4, False, "Strong", "BS4: Not evaluated", {})
    
    def _evaluate_BP1(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BP1, False, "Supporting", "BP1: Not evaluated", {})
    
    def _evaluate_BP2(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BP2, False, "Supporting", "BP2: Not evaluated", {})
    
    def _evaluate_BP3(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BP3, False, "Supporting", "BP3: Not evaluated", {})
    
    def _evaluate_BP4(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BP4, False, "Supporting", "BP4: Not evaluated", {})
    
    def _evaluate_BP5(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BP5, False, "Supporting", "BP5: Not evaluated", {})
    
    def _evaluate_BP6(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BP6, False, "Supporting", "BP6: Not evaluated", {})
    
    def _evaluate_BP7(self, variant: Variant) -> ACMGCriteriaResult:
        return ACMGCriteriaResult(ACMGCriteria.BP7, False, "Supporting", "BP7: Not evaluated", {})
    
    def _is_lof_variant(self, variant: Variant) -> bool:
        """Check if variant is predicted to cause loss of function."""
        # Check impact field
        if variant.impact:
            lof_terms = ['HIGH', 'stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant']
            return any(term in variant.impact.upper() for term in lof_terms)
        
        # Check HGVS notation for stop gain
        if variant.hgvs_p and 'Ter' in variant.hgvs_p:
            return True
        
        return False
    
    def _has_lof_mechanism(self, variant: Variant) -> bool:
        """Check if gene has known loss-of-function mechanism."""
        # This would typically query a database of genes with known LOF mechanisms
        # For now, return True for common tumor suppressor genes
        lof_genes = ['BRCA1', 'BRCA2', 'TP53', 'PTEN', 'APC', 'RB1', 'VHL']
        return variant.gene in lof_genes if variant.gene else False
    
    def _calculate_score(self, criteria_results: List[ACMGCriteriaResult]) -> int:
        """Calculate ACMG score from criteria results."""
        score = 0
        
        for result in criteria_results:
            if not result.met:
                continue
            
            if result.criteria in [ACMGCriteria.PVS1]:
                score += 8  # Very Strong
            elif result.criteria in [ACMGCriteria.PS1, ACMGCriteria.PS2, ACMGCriteria.PS3, ACMGCriteria.PS4]:
                score += 4  # Strong
            elif result.criteria in [ACMGCriteria.PM1, ACMGCriteria.PM2, ACMGCriteria.PM3, ACMGCriteria.PM4, ACMGCriteria.PM5, ACMGCriteria.PM6]:
                score += 2  # Moderate
            elif result.criteria in [ACMGCriteria.PP1, ACMGCriteria.PP2, ACMGCriteria.PP3, ACMGCriteria.PP4, ACMGCriteria.PP5]:
                score += 1  # Supporting
            elif result.criteria in [ACMGCriteria.BA1]:
                score -= 8  # Standalone Benign
            elif result.criteria in [ACMGCriteria.BS1, ACMGCriteria.BS2, ACMGCriteria.BS3, ACMGCriteria.BS4]:
                score -= 4  # Strong Benign
            elif result.criteria in [ACMGCriteria.BP1, ACMGCriteria.BP2, ACMGCriteria.BP3, ACMGCriteria.BP4, ACMGCriteria.BP5, ACMGCriteria.BP6, ACMGCriteria.BP7]:
                score -= 1  # Supporting Benign
        
        return score
    
    def _determine_classification(self, score: int) -> ACMGClassification:
        """Determine ACMG classification from score."""
        if score >= 8:
            return ACMGClassification.PATHOGENIC
        elif score >= 6:
            return ACMGClassification.LIKELY_PATHOGENIC
        elif score <= -8:
            return ACMGClassification.BENIGN
        elif score <= -6:
            return ACMGClassification.LIKELY_BENIGN
        else:
            return ACMGClassification.UNCERTAIN_SIGNIFICANCE
    
    def _generate_reasoning(self, criteria_results: List[ACMGCriteriaResult], classification: ACMGClassification) -> str:
        """Generate human-readable reasoning for classification."""
        met_criteria = [cr for cr in criteria_results if cr.met]
        
        if not met_criteria:
            return f"Classification: {classification.value}. No ACMG criteria met."
        
        reasoning_parts = [f"Classification: {classification.value}"]
        reasoning_parts.append("Applied criteria:")
        
        for result in met_criteria:
            reasoning_parts.append(f"  {result.criteria.value}: {result.reasoning}")
        
        return "\n".join(reasoning_parts)
    
    def _calculate_confidence(self, criteria_results: List[ACMGCriteriaResult]) -> float:
        """Calculate confidence score (0-1) based on criteria strength."""
        met_criteria = [cr for cr in criteria_results if cr.met]
        
        if not met_criteria:
            return 0.0
        
        # Calculate confidence based on criteria strength
        total_weight = 0
        weighted_sum = 0
        
        for result in met_criteria:
            if result.strength == "Very Strong":
                weight = 4
            elif result.strength == "Strong":
                weight = 3
            elif result.strength == "Moderate":
                weight = 2
            elif result.strength == "Supporting":
                weight = 1
            else:
                weight = 0
            
            total_weight += weight
            weighted_sum += weight
        
        return min(weighted_sum / max(total_weight, 1), 1.0) 