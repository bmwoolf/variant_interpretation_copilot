![Banner](assets/github_banner.png)

# Clinical Genomics Copilot (VCF → ACMG Report)

A modular command-line application for clinical genetic variant interpretation from VCF files. Built with Python for orchestration and Rust for high-performance VCF/BAM parsing.

## Architecture

- **Backend**: Python + Rust (PyO3 bindings)
- **VCF Parser**: Rust with `rust-htslib` for high-performance parsing
- **Annotation Engine**: Python modules for external API integration
- **ACMG Classifier**: Clinical variant classification system
- **Report Generator**: HTML/JSON output with rich formatting

## Quick Start

### Prerequisites
- Python 3.9+
- Rust 1.70+
- Cargo (Rust package manager)

### Installation

```bash
# Clone the repository
git clone https://github.com/bmwoolf/variant_interpretation_copilot
cd variant_interpretation_copilot

# Install Python dependencies
pip install -r requirements.txt

# Build Rust components
cd parser && cargo build --release
cd ..

# Install the CLI tool
pip install -e .
```

### Usage

```bash
# Basic usage
python -m vcf_copilot input.vcf --output report.html

# With phenotype filtering
python -m vcf_copilot input.vcf --output report.html --phenotype HP:0001250

# JSON output for programmatic use
python -m vcf_copilot input.vcf --output variants.json --format json
```

## Core Modules

### 1. VCF Parser (Rust)
- High-performance VCF parsing via `rust-htslib`
- GRCh38 coordinate system support
- Parallel variant processing

### 2. Annotation Engine (Python)
- **ClinVar**: Clinical significance and disease associations
- **gnomAD**: Population allele frequencies
- **Ensembl**: Gene and transcript information
- **Variant Impact**: VEP/SnpEff prediction parsing
- **In silico predictions**: CADD, PolyPhen, SIFT scores

### 3. ACMG Classifier
- ACMG criteria implementation (PVS1, PM2, PP3, etc.)
- Variant classification (Pathogenic, VUS, Benign)
- Evidence-based reasoning and scoring

### 4. Report Generator
- Rich HTML reports with clinical context
- JSON output for integration
- Phenotype-based gene filtering

## Example Output

```
Gene: BRCA1
Variant: c.68_69delAG
Classification: Pathogenic
ACMG Tags: PVS1, PM2
Evidence: Predicted loss-of-function in tumor suppressor gene
Population Frequency: 0.0001% (gnomAD)
Clinical Significance: Pathogenic (ClinVar)
```

## Running Tests
```bash
# Python tests
pytest tests/

# Rust tests
cd parser && cargo test
```

## Configuration

Environment variables for API access:
- `CLINVAR_API_KEY`: ClinVar API access
- `ENSEMBL_API_KEY`: Ensembl REST API access
- `CADD_API_KEY`: CADD API access

## Roadmap

- [ ] CNV support via ExomeDepth
- [ ] Mitochondrial variant detection
- [ ] Pharmacogenomics (PGx) flags
- [ ] Nextflow pipeline integration
- [ ] Cloud-native microservice deployment
- [ ] LLM-powered gene-disease mechanism summaries

## License

MIT.

## Project Structure
```
vcf_copilot/
├── main.py                 # CLI entry point
├── parser/                 # Rust VCF parser
├── annotators/            # External API integrations
├── scoring/               # ACMG classification
├── report/                # Output generation
├── tests/                 # Unit tests
└── data/                  # Test data
```