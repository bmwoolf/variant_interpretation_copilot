"""
gnomAD annotation module for population allele frequencies.
"""

import logging
import time
from typing import Dict, Optional
import requests

from ..models import Variant

logger = logging.getLogger(__name__)


class GnomADAnnotator:
    """gnomAD API annotator for population allele frequencies."""
    
    def __init__(self):
        """Initialize gnomAD annotator."""
        self.base_url = "https://gnomad.broadinstitute.org/api"
        self.logger = logging.getLogger(__name__)
        
        # Rate limiting
        self.last_request_time = 0
        self.min_request_interval = 0.2  # 200ms between requests
        
    def annotate(self, variant: Variant) -> Optional[Dict]:
        """
        Annotate variant with gnomAD data.
        
        Args:
            variant: Variant to annotate
            
        Returns:
            gnomAD annotation data or None if not found
        """
        try:
            # Rate limiting
            self._rate_limit()
            
            # Query gnomAD API
            variant_data = self._query_gnomad(variant)
            if not variant_data:
                return None
            
            return {
                'allele_frequency': variant_data.get('allele_frequency'),
                'allele_count': variant_data.get('allele_count'),
                'total_count': variant_data.get('total_count'),
                'homozygote_count': variant_data.get('homozygote_count'),
                'chromosome': variant_data.get('chromosome'),
                'position': variant_data.get('position'),
                'reference_genome': variant_data.get('reference_genome'),
            }
            
        except Exception as e:
            self.logger.warning(f"gnomAD annotation failed for {variant.chrom}:{variant.pos}: {e}")
            return None
    
    def _rate_limit(self) -> None:
        """Implement rate limiting for API requests."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        
        if time_since_last < self.min_request_interval:
            time.sleep(self.min_request_interval - time_since_last)
        
        self.last_request_time = time.time()
    
    def _query_gnomad(self, variant: Variant) -> Optional[Dict]:
        """
        Query gnomAD API for variant information.
        
        Args:
            variant: Variant to query
            
        Returns:
            gnomAD variant data or None if not found
        """
        try:
            # Build GraphQL query
            query = self._build_graphql_query(variant)
            
            # Make request
            response = requests.post(
                self.base_url,
                json={'query': query},
                headers={'Content-Type': 'application/json'},
                timeout=30
            )
            response.raise_for_status()
            
            data = response.json()
            
            # Extract variant data
            if 'data' in data and 'variant' in data['data']:
                return data['data']['variant']
            
            return None
            
        except Exception as e:
            self.logger.warning(f"gnomAD query failed: {e}")
            return None
    
    def _build_graphql_query(self, variant: Variant) -> str:
        """
        Build GraphQL query for gnomAD.
        
        Args:
            variant: Variant to build query for
            
        Returns:
            GraphQL query string
        """
        # Convert chromosome format if needed
        chrom = variant.chrom
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        
        query = f"""
        {{
            variant(chromosome: "{chrom}", position: {variant.pos}, reference_genome: GRCh38) {{
                chromosome
                position
                reference_genome
                reference_allele
                alternate_allele
                allele_frequency
                allele_count
                total_count
                homozygote_count
                populations {{
                    id
                    allele_frequency
                    allele_count
                    total_count
                }}
            }}
        }}
        """
        
        return query
    
    def get_population_frequencies(self, variant: Variant) -> Dict[str, float]:
        """
        Get population-specific allele frequencies.
        
        Args:
            variant: Variant to get frequencies for
            
        Returns:
            Dictionary of population frequencies
        """
        try:
            variant_data = self._query_gnomad(variant)
            if not variant_data or 'populations' not in variant_data:
                return {}
            
            frequencies = {}
            for population in variant_data['populations']:
                pop_id = population['id']
                af = population.get('allele_frequency')
                if af is not None:
                    frequencies[pop_id] = af
            
            return frequencies
            
        except Exception as e:
            self.logger.warning(f"Failed to get population frequencies: {e}")
            return {}
    
    def is_common_variant(self, variant: Variant, threshold: float = 0.01) -> bool:
        """
        Check if variant is common in population.
        
        Args:
            variant: Variant to check
            threshold: Allele frequency threshold (default: 1%)
            
        Returns:
            True if variant is common, False otherwise
        """
        try:
            variant_data = self._query_gnomad(variant)
            if not variant_data:
                return False
            
            af = variant_data.get('allele_frequency')
            return af is not None and af > threshold
            
        except Exception as e:
            self.logger.warning(f"Failed to check if variant is common: {e}")
            return False 