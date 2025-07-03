"""
Main annotation engine that coordinates all annotation sources.
"""

import logging
import requests
from typing import List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

from ..models import Variant
from .clinvar import ClinVarAnnotator
from .gnomad import GnomADAnnotator
from .ensembl import EnsemblAnnotator

logger = logging.getLogger(__name__)


class AnnotationEngine:
    """Main annotation engine that coordinates all annotation sources."""
    
    def __init__(self, max_workers: int = 4):
        """
        Initialize annotation engine.
        
        Args:
            max_workers: Maximum number of concurrent annotation workers
        """
        self.max_workers = max_workers
        self.logger = logging.getLogger(__name__)
        
        # Initialize annotators
        self.clinvar_annotator = ClinVarAnnotator()
        self.gnomad_annotator = GnomADAnnotator()
        self.ensembl_annotator = EnsemblAnnotator()
        
        self.logger.info("Annotation engine initialized")
    
    def annotate(self, variant: Variant) -> Variant:
        """
        Annotate a single variant with all available sources.
        
        Args:
            variant: Variant to annotate
            
        Returns:
            Annotated variant
        """
        try:
            # Create a copy of the variant for annotation
            annotated_variant = variant.copy(deep=True)
            
            # Annotate with ClinVar (with timeout and better error handling)
            try:
                clinvar_data = self.clinvar_annotator.annotate(variant)
                if clinvar_data:
                    annotated_variant.clinvar_significance = clinvar_data.get('clinical_significance')
                    annotated_variant.clinvar_diseases = clinvar_data.get('diseases', [])
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
                self.logger.warning(f"ClinVar API timeout/connection error for {variant.chrom}:{variant.pos}: {e}")
            except Exception as e:
                self.logger.warning(f"ClinVar annotation failed for {variant.chrom}:{variant.pos}: {e}")
            
            # Annotate with gnomAD (with timeout and better error handling)
            try:
                gnomad_data = self.gnomad_annotator.annotate(variant)
                if gnomad_data:
                    annotated_variant.gnomad_af = gnomad_data.get('allele_frequency')
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
                self.logger.warning(f"gnomAD API timeout/connection error for {variant.chrom}:{variant.pos}: {e}")
            except Exception as e:
                self.logger.warning(f"gnomAD annotation failed for {variant.chrom}:{variant.pos}: {e}")
            
            # Annotate with Ensembl (with timeout and better error handling)
            try:
                ensembl_data = self.ensembl_annotator.annotate(variant)
                if ensembl_data:
                    if not annotated_variant.gene:
                        annotated_variant.gene = ensembl_data.get('gene_name')
                    if not annotated_variant.transcript:
                        annotated_variant.transcript = ensembl_data.get('transcript_id')
                    if not annotated_variant.hgvs_c:
                        annotated_variant.hgvs_c = ensembl_data.get('hgvs_c')
                    if not annotated_variant.hgvs_p:
                        annotated_variant.hgvs_p = ensembl_data.get('hgvs_p')
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
                self.logger.warning(f"Ensembl API timeout/connection error for {variant.chrom}:{variant.pos}: {e}")
            except Exception as e:
                self.logger.warning(f"Ensembl annotation failed for {variant.chrom}:{variant.pos}: {e}")
            
            # Extract in silico predictions from INFO field
            self._extract_in_silico_predictions(annotated_variant)
            
            return annotated_variant
            
        except Exception as e:
            self.logger.error(f"Annotation failed for variant {variant.chrom}:{variant.pos}: {e}")
            return variant
    
    def annotate_batch(self, variants: List[Variant]) -> List[Variant]:
        """
        Annotate multiple variants in parallel.
        
        Args:
            variants: List of variants to annotate
            
        Returns:
            List of annotated variants
        """
        annotated_variants = []
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit annotation tasks
            future_to_variant = {
                executor.submit(self.annotate, variant): variant 
                for variant in variants
            }
            
            # Collect results
            for future in as_completed(future_to_variant):
                try:
                    annotated_variant = future.result()
                    annotated_variants.append(annotated_variant)
                except Exception as e:
                    variant = future_to_variant[future]
                    self.logger.error(f"Batch annotation failed for {variant.chrom}:{variant.pos}: {e}")
                    annotated_variants.append(variant)  # Return original variant
        
        return annotated_variants
    
    def _extract_in_silico_predictions(self, variant: Variant) -> None:
        """
        Extract in silico prediction scores from VCF INFO field.
        
        Args:
            variant: Variant to extract predictions from
        """
        info = variant.info
        
        # CADD score
        if 'CADD_RAW' in info:
            try:
                variant.cadd_score = float(info['CADD_RAW'])
            except (ValueError, TypeError):
                pass
        
        # PolyPhen score
        if 'PolyPhen' in info:
            try:
                polyphen_str = str(info['PolyPhen'])
                if ':' in polyphen_str:
                    score_str = polyphen_str.split(':')[1]
                    variant.polyphen_score = float(score_str)
            except (ValueError, TypeError):
                pass
        
        # SIFT score
        if 'SIFT' in info:
            try:
                sift_str = str(info['SIFT'])
                if ':' in sift_str:
                    score_str = sift_str.split(':')[1]
                    variant.sift_score = float(score_str)
            except (ValueError, TypeError):
                pass
    
    def get_annotation_summary(self, variants: List[Variant]) -> dict:
        """
        Get summary statistics for annotation coverage.
        
        Args:
            variants: List of annotated variants
            
        Returns:
            Dictionary with annotation statistics
        """
        total_variants = len(variants)
        
        stats = {
            'total_variants': total_variants,
            'clinvar_annotated': sum(1 for v in variants if v.clinvar_significance),
            'gnomad_annotated': sum(1 for v in variants if v.gnomad_af is not None),
            'ensembl_annotated': sum(1 for v in variants if v.gene),
            'cadd_annotated': sum(1 for v in variants if v.cadd_score is not None),
            'polyphen_annotated': sum(1 for v in variants if v.polyphen_score is not None),
            'sift_annotated': sum(1 for v in variants if v.sift_score is not None),
        }
        
        # Calculate percentages
        for key in stats:
            if key != 'total_variants' and total_variants > 0:
                stats[f'{key}_percentage'] = (stats[key] / total_variants) * 100
        
        return stats 