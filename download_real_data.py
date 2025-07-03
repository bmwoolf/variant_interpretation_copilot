#!/usr/bin/env python3
"""
Download real VCF data for testing the variant interpretation pipeline.
"""

import os
import requests
import gzip
import shutil
from pathlib import Path
from urllib.parse import urlparse

def download_file(url: str, output_path: Path, chunk_size: int = 8192) -> bool:
    """Download a file with progress indication."""
    try:
        print(f"Downloading {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"\rProgress: {percent:.1f}%", end='', flush=True)
        
        print(f"\nDownloaded: {output_path}")
        return True
        
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return False

def extract_gzip(gzip_path: Path, output_path: Path) -> bool:
    """Extract a gzipped file."""
    try:
        print(f"Extracting {gzip_path}...")
        with gzip.open(gzip_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Extracted: {output_path}")
        return True
    except Exception as e:
        print(f"Error extracting {gzip_path}: {e}")
        return False

def create_data_directory():
    """Create data directory if it doesn't exist."""
    data_dir = Path("data/real")
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir

def download_clinvar_variants():
    """Download a small subset of ClinVar variants."""
    data_dir = create_data_directory()
    
    # ClinVar pathogenic variants (small subset)
    clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    clinvar_gz = data_dir / "clinvar.vcf.gz"
    clinvar_vcf = data_dir / "clinvar_subset.vcf"
    
    if download_file(clinvar_url, clinvar_gz):
        # Extract and create a small subset
        if extract_gzip(clinvar_gz, clinvar_vcf):
            # Create a smaller subset (first 50 variants)
            subset_path = data_dir / "clinvar_small.vcf"
            header_lines = []
            variant_lines = []
            variant_count = 0
            max_variants = 50
            
            with open(clinvar_vcf, 'r') as f_in:
                for line in f_in:
                    if line.startswith('#'):
                        header_lines.append(line)
                    elif variant_count < max_variants:
                        variant_lines.append(line)
                        variant_count += 1
                    else:
                        break
            
            # Write the subset file
            with open(subset_path, 'w') as f_out:
                f_out.writelines(header_lines)
                f_out.writelines(variant_lines)
            
            print(f"Created ClinVar subset: {subset_path} ({variant_count} variants)")
            
            # Clean up large files
            clinvar_gz.unlink(missing_ok=True)
            clinvar_vcf.unlink(missing_ok=True)

def download_1000genomes_variants():
    """Download a small subset of 1000 Genomes variants."""
    data_dir = create_data_directory()
    
    # 1000 Genomes phase 3 variants (chromosome 22 as example)
    kg_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    kg_gz = data_dir / "1000genomes_chr22.vcf.gz"
    kg_vcf = data_dir / "1000genomes_chr22.vcf"
    
    if download_file(kg_url, kg_gz):
        if extract_gzip(kg_gz, kg_vcf):
            # Create a smaller subset
            subset_path = data_dir / "1000genomes_small.vcf"
            header_lines = []
            variant_lines = []
            variant_count = 0
            max_variants = 100
            
            with open(kg_vcf, 'r') as f_in:
                for line in f_in:
                    if line.startswith('#'):
                        header_lines.append(line)
                    elif variant_count < max_variants:
                        variant_lines.append(line)
                        variant_count += 1
                    else:
                        break
            
            # Write the subset file
            with open(subset_path, 'w') as f_out:
                f_out.writelines(header_lines)
                f_out.writelines(variant_lines)
            
            print(f"Created 1000 Genomes subset: {subset_path} ({variant_count} variants)")
            
            # Clean up large files
            kg_gz.unlink(missing_ok=True)
            kg_vcf.unlink(missing_ok=True)

def download_gnomad_variants():
    """Download a small subset of gnomAD variants."""
    data_dir = create_data_directory()
    
    # gnomAD exomes (chromosome 1 as example)
    gnomad_url = "https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.chr1.vcf.bgz"
    gnomad_bgz = data_dir / "gnomad_chr1.vcf.bgz"
    gnomad_vcf = data_dir / "gnomad_chr1.vcf"
    
    if download_file(gnomad_url, gnomad_bgz):
        # Note: .bgz files need special handling, but for now we'll skip this
        # as it requires additional tools like tabix
        print("gnomAD .bgz files require tabix for processing. Skipping for now.")
        gnomad_bgz.unlink(missing_ok=True)

def main():
    """Download real VCF data for testing."""
    print("🧬 Downloading real VCF data for testing...")
    print("=" * 50)
    
    # Create data directory
    data_dir = create_data_directory()
    print(f"Data will be saved to: {data_dir}")
    print()
    
    # Download different datasets
    print("1. Downloading ClinVar variants...")
    download_clinvar_variants()
    print()
    
    print("2. Downloading 1000 Genomes variants...")
    download_1000genomes_variants()
    print()
    
    print("3. Downloading gnomAD variants...")
    download_gnomad_variants()
    print()
    
    print("✅ Download complete!")
    print("\nAvailable test files:")
    for file in data_dir.glob("*.vcf"):
        print(f"  - {file}")
    
    print("\nYou can now test with real data:")
    print(f"python -m vcf_copilot main {data_dir}/clinvar_small.vcf --output real_report.html")

if __name__ == "__main__":
    main() 