"""
ClinVar annotation module for clinical significance and disease associations.
"""

import logging
import time
from typing import Dict, Optional, List
import requests

from ..models import Variant

logger = logging.getLogger(__name__)


class ClinVarAnnotator:
    """ClinVar REST API annotator for clinical significance."""
    
    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize ClinVar annotator.
        
        Args:
            api_key: ClinVar API key (optional)
        """
        self.api_key = api_key
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.logger = logging.getLogger(__name__)
        
        # Rate limiting
        self.last_request_time = 0
        self.min_request_interval = 0.1  # 100ms between requests
        
    def annotate(self, variant: Variant) -> Optional[Dict]:
        """
        Annotate variant with ClinVar data.
        
        Args:
            variant: Variant to annotate
            
        Returns:
            ClinVar annotation data or None if not found
        """
        try:
            # Rate limiting
            self._rate_limit()
            
            # Search for variant in ClinVar
            variant_id = self._search_variant(variant)
            if not variant_id:
                return None
            
            # Get detailed variant information
            variant_data = self._get_variant_details(variant_id)
            if not variant_data:
                return None
            
            return {
                'clinical_significance': variant_data.get('clinical_significance'),
                'diseases': variant_data.get('diseases', []),
                'review_status': variant_data.get('review_status'),
                'last_evaluated': variant_data.get('last_evaluated'),
                'variant_id': variant_id,
            }
            
        except Exception as e:
            self.logger.warning(f"ClinVar annotation failed for {variant.chrom}:{variant.pos}: {e}")
            return None
    
    def _rate_limit(self) -> None:
        """Implement rate limiting for API requests."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        
        if time_since_last < self.min_request_interval:
            time.sleep(self.min_request_interval - time_since_last)
        
        self.last_request_time = time.time()
    
    def _search_variant(self, variant: Variant) -> Optional[str]:
        """
        Search for variant in ClinVar database.
        
        Args:
            variant: Variant to search for
            
        Returns:
            ClinVar variant ID or None if not found
        """
        try:
            # Build search query
            query = self._build_search_query(variant)
            
            # Search parameters
            params = {
                'db': 'clinvar',
                'term': query,
                'retmode': 'json',
                'retmax': 1,
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            # Make request
            response = requests.get(f"{self.base_url}/esearch.fcgi", params=params)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract variant ID
            if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                id_list = data['esearchresult']['idlist']
                if id_list:
                    return id_list[0]
            
            return None
            
        except Exception as e:
            self.logger.warning(f"ClinVar search failed: {e}")
            return None
    
    def _build_search_query(self, variant: Variant) -> str:
        """
        Build search query for ClinVar.
        
        Args:
            variant: Variant to build query for
            
        Returns:
            Search query string
        """
        # Basic query with chromosome and position
        query_parts = [
            f"{variant.chrom}[chr]",
            f"{variant.pos}[pos]",
        ]
        
        # Add gene if available
        if variant.gene:
            query_parts.append(f"{variant.gene}[gene]")
        
        # Add HGVS notation if available
        if variant.hgvs_c:
            query_parts.append(f'"{variant.hgvs_c}"[hgvs]')
        
        if variant.hgvs_p:
            query_parts.append(f'"{variant.hgvs_p}"[hgvs]')
        
        return " AND ".join(query_parts)
    
    def _get_variant_details(self, variant_id: str) -> Optional[Dict]:
        """
        Get detailed variant information from ClinVar.
        
        Args:
            variant_id: ClinVar variant ID
            
        Returns:
            Variant details or None if not found
        """
        try:
            # Get summary
            params = {
                'db': 'clinvar',
                'id': variant_id,
                'retmode': 'json',
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(f"{self.base_url}/esummary.fcgi", params=params)
            response.raise_for_status()
            
            data = response.json()
            
            if 'result' not in data or variant_id not in data['result']:
                return None
            
            variant_data = data['result'][variant_id]
            
            # Extract relevant information
            return {
                'clinical_significance': variant_data.get('clinical_significance'),
                'diseases': self._extract_diseases(variant_data),
                'review_status': variant_data.get('review_status'),
                'last_evaluated': variant_data.get('last_evaluated'),
            }
            
        except Exception as e:
            self.logger.warning(f"Failed to get ClinVar details for {variant_id}: {e}")
            return None
    
    def _extract_diseases(self, variant_data: Dict) -> List[str]:
        """
        Extract disease names from ClinVar variant data.
        
        Args:
            variant_data: ClinVar variant data
            
        Returns:
            List of disease names
        """
        diseases = []
        
        try:
            # Extract from phenotype information
            if 'phenotype_ids' in variant_data:
                for phenotype in variant_data['phenotype_ids']:
                    if 'disease_name' in phenotype:
                        diseases.append(phenotype['disease_name'])
            
            # Extract from trait information
            if 'trait_set' in variant_data:
                for trait in variant_data['trait_set']:
                    if 'trait_name' in trait:
                        diseases.append(trait['trait_name'])
            
        except Exception as e:
            self.logger.warning(f"Failed to extract diseases: {e}")
        
        return list(set(diseases))  # Remove duplicates 