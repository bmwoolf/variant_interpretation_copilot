"""
Ensembl annotation module for gene and transcript information.
"""

import logging
import time
from typing import Dict, Optional
import requests

from ..models import Variant

logger = logging.getLogger(__name__)


class EnsemblAnnotator:
    """Ensembl REST API annotator for gene and transcript information."""
    
    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize Ensembl annotator.
        
        Args:
            api_key: Ensembl API key (optional)
        """
        self.api_key = api_key
        self.base_url = "https://rest.ensembl.org"
        self.logger = logging.getLogger(__name__)
        
        # Rate limiting
        self.last_request_time = 0
        self.min_request_interval = 0.1  # 100ms between requests
        
    def annotate(self, variant: Variant) -> Optional[Dict]:
        """
        Annotate variant with Ensembl data.
        
        Args:
            variant: Variant to annotate
            
        Returns:
            Ensembl annotation data or None if not found
        """
        try:
            # Rate limiting
            self._rate_limit()
            
            # Get variant consequences
            consequences = self._get_variant_consequences(variant)
            if not consequences:
                return None
            
            # Extract relevant information
            return self._extract_annotation_data(consequences)
            
        except Exception as e:
            self.logger.warning(f"Ensembl annotation failed for {variant.chrom}:{variant.pos}: {e}")
            return None
    
    def _rate_limit(self) -> None:
        """Implement rate limiting for API requests."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        
        if time_since_last < self.min_request_interval:
            time.sleep(self.min_request_interval - time_since_last)
        
        self.last_request_time = time.time()
    
    def _get_variant_consequences(self, variant: Variant) -> Optional[Dict]:
        """
        Get variant consequences from Ensembl VEP API.
        
        Args:
            variant: Variant to get consequences for
            
        Returns:
            Variant consequences data or None if not found
        """
        try:
            # Build variant identifier
            variant_id = self._build_variant_id(variant)
            
            # Query VEP API
            url = f"{self.base_url}/vep/human/id/{variant_id}"
            
            headers = {
                'Content-Type': 'application/json',
            }
            
            if self.api_key:
                headers['Authorization'] = f'Bearer {self.api_key}'
            
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            if data and len(data) > 0:
                return data[0]  # Return first result
            
            return None
            
        except Exception as e:
            self.logger.warning(f"Ensembl VEP query failed: {e}")
            return None
    
    def _build_variant_id(self, variant: Variant) -> str:
        """
        Build Ensembl variant identifier.
        
        Args:
            variant: Variant to build ID for
            
        Returns:
            Ensembl variant ID
        """
        # Convert chromosome format
        chrom = variant.chrom
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        
        # Build variant ID in format: chr:pos:ref:alt
        return f"{chrom}:{variant.pos}:{variant.ref}:{variant.alt}"
    
    def _extract_annotation_data(self, consequences: Dict) -> Dict:
        """
        Extract relevant annotation data from consequences.
        
        Args:
            consequences: Variant consequences data
            
        Returns:
            Extracted annotation data
        """
        annotation_data = {}
        
        try:
            # Get transcript consequences
            transcript_consequences = consequences.get('transcript_consequences', [])
            
            if transcript_consequences:
                # Use first transcript consequence
                transcript = transcript_consequences[0]
                
                annotation_data.update({
                    'gene_name': transcript.get('gene_symbol'),
                    'transcript_id': transcript.get('transcript_id'),
                    'gene_id': transcript.get('gene_id'),
                    'consequence_terms': transcript.get('consequence_terms', []),
                    'impact': transcript.get('impact'),
                    'biotype': transcript.get('biotype'),
                })
                
                # Get HGVS notations
                if 'hgvs' in transcript:
                    annotation_data['hgvs_c'] = transcript['hgvs']
                
                if 'hgvsp' in transcript:
                    annotation_data['hgvs_p'] = transcript['hgvsp']
                
                # Get protein information
                if 'protein_start' in transcript:
                    annotation_data['protein_start'] = transcript['protein_start']
                
                if 'protein_end' in transcript:
                    annotation_data['protein_end'] = transcript['protein_end']
                
                # Get cDNA information
                if 'cdna_start' in transcript:
                    annotation_data['cdna_start'] = transcript['cdna_start']
                
                if 'cdna_end' in transcript:
                    annotation_data['cdna_end'] = transcript['cdna_end']
            
            # Get regulatory consequences
            regulatory_consequences = consequences.get('regulatory_consequences', [])
            if regulatory_consequences:
                annotation_data['regulatory_consequences'] = regulatory_consequences
            
        except Exception as e:
            self.logger.warning(f"Failed to extract annotation data: {e}")
        
        return annotation_data
    
    def get_gene_info(self, gene_id: str) -> Optional[Dict]:
        """
        Get detailed gene information from Ensembl.
        
        Args:
            gene_id: Ensembl gene ID
            
        Returns:
            Gene information or None if not found
        """
        try:
            url = f"{self.base_url}/lookup/id/{gene_id}"
            
            headers = {
                'Content-Type': 'application/json',
            }
            
            if self.api_key:
                headers['Authorization'] = f'Bearer {self.api_key}'
            
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()
            
            return response.json()
            
        except Exception as e:
            self.logger.warning(f"Failed to get gene info for {gene_id}: {e}")
            return None
    
    def get_transcript_info(self, transcript_id: str) -> Optional[Dict]:
        """
        Get detailed transcript information from Ensembl.
        
        Args:
            transcript_id: Ensembl transcript ID
            
        Returns:
            Transcript information or None if not found
        """
        try:
            url = f"{self.base_url}/lookup/id/{transcript_id}"
            
            headers = {
                'Content-Type': 'application/json',
            }
            
            if self.api_key:
                headers['Authorization'] = f'Bearer {self.api_key}'
            
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()
            
            return response.json()
            
        except Exception as e:
            self.logger.warning(f"Failed to get transcript info for {transcript_id}: {e}")
            return None 