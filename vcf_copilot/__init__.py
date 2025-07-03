"""
Clinical Genomics Copilot for VCF variant interpretation.

A modular command-line application for clinical genetic variant interpretation
from VCF files using Python for orchestration and Rust for high-performance parsing.
"""

__version__ = "0.1.0"
__author__ = "Clinical Genomics Team"

from .models import Variant, Annotation, ACMGClassification
from .scoring import ACMGClassifier
from .report import ReportGenerator

__all__ = [
    "Variant",
    "Annotation", 
    "ACMGClassification",
    "ACMGClassifier",
    "ReportGenerator",
] 