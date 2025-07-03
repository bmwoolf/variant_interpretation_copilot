"""
Core data models for the VCF Copilot application.
"""

from enum import Enum
from typing import Dict, List, Optional, Any
from pydantic import BaseModel, Field


class VariantType(str, Enum):
    """Variant types."""
    SNV = "SNV"
    INDEL = "INDEL"
    CNV = "CNV"
    SV = "SV"


class ACMGClassification(str, Enum):
    """ACMG variant classifications."""
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    UNCERTAIN_SIGNIFICANCE = "Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


class ACMGCriteria(str, Enum):
    """ACMG criteria codes."""
    # Pathogenic Very Strong
    PVS1 = "PVS1"
    # Pathogenic Strong
    PS1 = "PS1"
    PS2 = "PS2"
    PS3 = "PS3"
    PS4 = "PS4"
    # Pathogenic Moderate
    PM1 = "PM1"
    PM2 = "PM2"
    PM3 = "PM3"
    PM4 = "PM4"
    PM5 = "PM5"
    PM6 = "PM6"
    # Pathogenic Supporting
    PP1 = "PP1"
    PP2 = "PP2"
    PP3 = "PP3"
    PP4 = "PP4"
    PP5 = "PP5"
    # Benign Standalone
    BA1 = "BA1"
    # Benign Strong
    BS1 = "BS1"
    BS2 = "BS2"
    BS3 = "BS3"
    BS4 = "BS4"
    # Benign Supporting
    BP1 = "BP1"
    BP2 = "BP2"
    BP3 = "BP3"
    BP4 = "BP4"
    BP5 = "BP5"
    BP6 = "BP6"
    BP7 = "BP7"


class Variant(BaseModel):
    """Core variant model."""
    chrom: str = Field(..., description="Chromosome")
    pos: int = Field(..., description="Position (1-based)")
    ref: str = Field(..., description="Reference allele")
    alt: str = Field(..., description="Alternate allele")
    variant_type: VariantType = Field(..., description="Type of variant")
    
    # Optional fields from VCF
    id: Optional[str] = Field(None, description="Variant ID")
    qual: Optional[float] = Field(None, description="Quality score")
    filter: Optional[str] = Field(None, description="Filter status")
    info: Dict[str, Any] = Field(default_factory=dict, description="INFO field")
    
    # Annotations
    gene: Optional[str] = Field(None, description="Gene symbol")
    transcript: Optional[str] = Field(None, description="Transcript ID")
    hgvs_c: Optional[str] = Field(None, description="cDNA notation")
    hgvs_p: Optional[str] = Field(None, description="Protein notation")
    impact: Optional[str] = Field(None, description="Variant impact")
    
    # Clinical annotations
    clinvar_significance: Optional[str] = Field(None, description="ClinVar clinical significance")
    clinvar_diseases: List[str] = Field(default_factory=list, description="ClinVar disease associations")
    gnomad_af: Optional[float] = Field(None, description="gnomAD allele frequency")
    cadd_score: Optional[float] = Field(None, description="CADD score")
    polyphen_score: Optional[float] = Field(None, description="PolyPhen score")
    sift_score: Optional[float] = Field(None, description="SIFT score")
    
    # ACMG classification
    acmg_classification: Optional[ACMGClassification] = Field(None, description="ACMG classification")
    acmg_criteria: List[ACMGCriteria] = Field(default_factory=list, description="Applied ACMG criteria")
    acmg_score: Optional[int] = Field(None, description="ACMG score")
    acmg_reasoning: Optional[str] = Field(None, description="ACMG reasoning")
    
    class Config:
        use_enum_values = True


class Annotation(BaseModel):
    """Annotation result model."""
    variant_id: str = Field(..., description="Variant identifier")
    source: str = Field(..., description="Annotation source")
    data: Dict[str, Any] = Field(..., description="Annotation data")
    timestamp: Optional[str] = Field(None, description="Annotation timestamp")


class ACMGResult(BaseModel):
    """ACMG classification result."""
    variant: Variant = Field(..., description="Variant")
    classification: ACMGClassification = Field(..., description="ACMG classification")
    criteria: List[ACMGCriteria] = Field(..., description="Applied criteria")
    score: int = Field(..., description="ACMG score")
    reasoning: str = Field(..., description="Detailed reasoning")
    confidence: float = Field(..., description="Confidence score (0-1)")


class ReportConfig(BaseModel):
    """Report generation configuration."""
    output_format: str = Field("html", description="Output format (html/json)")
    include_evidence: bool = Field(True, description="Include evidence details")
    include_recommendations: bool = Field(True, description="Include clinical recommendations")
    phenotype_filter: Optional[str] = Field(None, description="HPO phenotype filter")
    min_quality: Optional[float] = Field(None, description="Minimum quality score")
    max_gnomad_af: Optional[float] = Field(0.01, description="Maximum gnomAD AF for reporting") 