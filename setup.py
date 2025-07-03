#!/usr/bin/env python3
"""
Setup script for VCF Copilot.
"""

from setuptools import setup, find_packages

if __name__ == "__main__":
    setup(
        name="vcf-copilot",
        version="0.1.0",
        description="Clinical Genomics Copilot for VCF variant interpretation",
        author="Clinical Genomics Team",
        author_email="team@example.com",
        packages=find_packages(),
        install_requires=[
            "typer>=0.9.0",
            "rich>=13.0.0",
            "pydantic>=2.0.0",
            "jinja2>=3.1.0",
            "requests>=2.31.0",
            "httpx>=0.24.0",
            "cyvcf2>=0.30.0",
            "pysam>=0.21.0",
            "pandas>=2.0.0",
            "numpy>=1.24.0",
        ],
        extras_require={
            "dev": [
                "pytest>=7.4.0",
                "pytest-cov>=4.1.0",
                "black>=23.0.0",
                "isort>=5.12.0",
                "flake8>=6.0.0",
            ],
            "llm": [
                "openai>=1.0.0",
            ],
        },
        entry_points={
            "console_scripts": [
                "vcf-copilot=vcf_copilot.main:app",
            ],
        },
        python_requires=">=3.9",
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    ) 