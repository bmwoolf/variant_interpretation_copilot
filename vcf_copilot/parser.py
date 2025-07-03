"""
VCF parser module for high-performance variant parsing.
"""

import logging
from pathlib import Path
from typing import List, Optional, Dict, Any
from dataclasses import dataclass

import cyvcf2
from rich.console import Console

from .models import Variant, VariantType

logger = logging.getLogger(__name__)
console = Console()


@dataclass
class ValidationResult:
    """VCF validation result."""
    is_valid: bool
    variant_count: int
    sample_count: int
    errors: List[str]


class VCFParser:
    """High-performance VCF parser using cyvcf2."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def parse(self, vcf_path: Path) -> List[Variant]:
        """
        Parse VCF file and return list of Variant objects.
        
        Args:
            vcf_path: Path to VCF file
            
        Returns:
            List of parsed Variant objects
        """
        variants = []
        
        try:
            vcf_reader = cyvcf2.VCF(str(vcf_path))
            
            for record in vcf_reader:
                try:
                    variant = self._parse_record(record)
                    if variant:
                        variants.append(variant)
                except Exception as e:
                    self.logger.warning(f"Skipping malformed record at {getattr(record, 'CHROM', '?')}:{getattr(record, 'POS', '?')}: {e}")
                    
        except Exception as e:
            self.logger.error(f"Error parsing VCF file: {e}")
            raise
            
        self.logger.info(f"Parsed {len(variants)} variants from {vcf_path}")
        return variants
    
    def _parse_record(self, record) -> Optional[Variant]:
        """
        Parse a single VCF record into a Variant object.
        
        Args:
            record: cyvcf2 record object
            
        Returns:
            Variant object or None if parsing fails
        """
        try:
            # Determine variant type
            variant_type = self._determine_variant_type(record)
            
            # Parse INFO field with better error handling
            try:
                info = self._parse_info(record)
            except Exception as e:
                self.logger.warning(f"Failed to parse INFO for record at {getattr(record, 'CHROM', '?')}:{getattr(record, 'POS', '?')}: {e}")
                info = {}  # Use empty dict as fallback
            
            # Extract gene and transcript information
            gene = info.get('Gene_Name', info.get('Gene', None))
            transcript = info.get('Transcript_ID', info.get('Feature', None))
            hgvs_c = info.get('HGVSc', None)
            hgvs_p = info.get('HGVSp', None)
            impact = info.get('IMPACT', info.get('Consequence', None))
            
            variant = Variant(
                chrom=record.CHROM,
                pos=record.POS,
                ref=record.REF,
                alt=str(record.ALT[0]) if record.ALT else "",
                variant_type=variant_type,
                id=record.ID if record.ID else None,
                qual=record.QUAL if record.QUAL else None,
                filter=record.FILTER if record.FILTER else None,
                info=info,
                gene=gene,
                transcript=transcript,
                hgvs_c=hgvs_c,
                hgvs_p=hgvs_p,
                impact=impact,
            )
            
            return variant
            
        except Exception as e:
            self.logger.warning(f"Failed to parse variant at {record.CHROM}:{record.POS}: {e}")
            return None
    
    def _determine_variant_type(self, record) -> VariantType:
        """Determine variant type from VCF record."""
        ref_len = len(record.REF)
        alt_len = len(str(record.ALT[0])) if record.ALT else 0
        
        if ref_len == 1 and alt_len == 1:
            return VariantType.SNV
        elif ref_len != alt_len:
            return VariantType.INDEL
        else:
            # Check for structural variants in INFO
            sv_type = record.INFO.get('SVTYPE')
            if sv_type:
                return VariantType.SV
            else:
                return VariantType.SNV  # Default fallback
    
    def _parse_info(self, record) -> Dict[str, Any]:
        """Parse INFO field into dictionary."""
        info = {}
        
        try:
            # Try to access INFO field safely
            if hasattr(record, 'INFO') and record.INFO:
                # Use dict() to convert to a regular Python dict
                raw_info = dict(record.INFO)
                
                for key, value in raw_info.items():
                    # Handle different value types properly
                    if isinstance(value, (list, tuple)):
                        # Convert each element to string and join
                        try:
                            info[key] = ','.join(str(v) for v in value)
                        except Exception as e:
                            # If individual elements can't be converted, use repr
                            self.logger.debug(f"Failed to convert tuple/list element: {e}")
                            info[key] = repr(value)
                    elif isinstance(value, bytes):
                        # Handle bytes objects
                        info[key] = value.decode('utf-8', errors='ignore')
                    elif hasattr(value, 'encode'):
                        # Handle objects with encode method (like bytes-like objects)
                        try:
                            info[key] = str(value)
                        except Exception:
                            info[key] = repr(value)
                    else:
                        info[key] = value
        except Exception as e:
            self.logger.warning(f"Failed to parse INFO for record at {getattr(record, 'CHROM', '?')}:{getattr(record, 'POS', '?')}: {e}")
        
        return info
    
    def validate(self, vcf_path: Path) -> ValidationResult:
        """
        Validate VCF file format and content.
        
        Args:
            vcf_path: Path to VCF file
            
        Returns:
            ValidationResult object
        """
        errors = []
        variant_count = 0
        sample_count = 0
        
        try:
            vcf_reader = cyvcf2.VCF(str(vcf_path))
            
            # Check header
            if not vcf_reader.samples:
                errors.append("No samples found in VCF header")
            else:
                sample_count = len(vcf_reader.samples)
            
            # Check for required header fields
            required_headers = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
            header_fields = [field for field in vcf_reader.raw_header.split('\n') if field.startswith('#CHROM')]
            
            if not header_fields:
                errors.append("Missing #CHROM header line")
            
            # Count variants and check for parsing errors
            for record in vcf_reader:
                variant_count += 1
                
                # Basic validation checks
                if not record.CHROM:
                    errors.append(f"Missing chromosome at variant {variant_count}")
                
                if not record.POS or record.POS <= 0:
                    errors.append(f"Invalid position at variant {variant_count}")
                
                if not record.REF:
                    errors.append(f"Missing reference allele at variant {variant_count}")
                
                if not record.ALT:
                    errors.append(f"Missing alternate allele at variant {variant_count}")
                
                # Limit error collection to avoid overwhelming output
                if len(errors) >= 10:
                    errors.append("... (additional errors truncated)")
                    break
                    
        except Exception as e:
            errors.append(f"Failed to read VCF file: {e}")
        
        is_valid = len(errors) == 0
        
        return ValidationResult(
            is_valid=is_valid,
            variant_count=variant_count,
            sample_count=sample_count,
            errors=errors
        )
    
    def get_header_info(self, vcf_path: Path) -> Dict[str, Any]:
        """
        Extract VCF header information.
        
        Args:
            vcf_path: Path to VCF file
            
        Returns:
            Dictionary with header information
        """
        try:
            vcf_reader = cyvcf2.VCF(str(vcf_path))
            
            header_info = {
                'samples': vcf_reader.samples,
                'sample_count': len(vcf_reader.samples),
                'raw_header': vcf_reader.raw_header,
            }
            
            return header_info
            
        except Exception as e:
            self.logger.error(f"Error reading VCF header: {e}")
            raise 