#!/usr/bin/env python3
"""
Example usage of VCF Copilot for clinical variant interpretation.
"""

import sys
from pathlib import Path

# Add the current directory to Python path for imports
sys.path.insert(0, str(Path(__file__).parent))

from vcf_copilot.parser import VCFParser
from vcf_copilot.annotators.engine import AnnotationEngine
from vcf_copilot.scoring import ACMGClassifier
from vcf_copilot.report import ReportGenerator
from vcf_copilot.models import ReportConfig


def main():
    """Example workflow for variant interpretation."""
    print("🧬 VCF Copilot - Clinical Variant Interpretation Example")
    print("=" * 60)
    
    # Check if test VCF exists
    vcf_path = Path("data/test.vcf")
    if not vcf_path.exists():
        print(f"❌ Test VCF file not found: {vcf_path}")
        print("Please ensure the test VCF file exists in the data/ directory.")
        return
    
    try:
        # Step 1: Parse VCF file
        print("\n📖 Step 1: Parsing VCF file...")
        parser = VCFParser()
        variants = parser.parse(vcf_path)
        print(f"✅ Parsed {len(variants)} variants")
        
        # Step 2: Annotate variants
        print("\n🔍 Step 2: Annotating variants...")
        annotator = AnnotationEngine()
        annotated_variants = []
        
        for i, variant in enumerate(variants, 1):
            print(f"  Annotating variant {i}/{len(variants)}: {variant.chrom}:{variant.pos}")
            annotated_variant = annotator.annotate(variant)
            annotated_variants.append(annotated_variant)
        
        print("✅ Annotation complete")
        
        # Step 3: Classify variants
        print("\n🎯 Step 3: Classifying variants...")
        classifier = ACMGClassifier()
        classified_variants = []
        
        for i, variant in enumerate(annotated_variants, 1):
            print(f"  Classifying variant {i}/{len(annotated_variants)}: {variant.chrom}:{variant.pos}")
            result = classifier.classify(variant)
            classified_variants.append(result)
        
        print("✅ Classification complete")
        
        # Step 4: Generate report
        print("\n📊 Step 4: Generating report...")
        report_generator = ReportGenerator()
        config = ReportConfig(
            output_format="html",
            include_evidence=True,
            include_recommendations=True,
            max_gnomad_af=0.01
        )
        
        output_path = Path("example_report.html")
        report_generator.generate(classified_variants, output_path, config)
        print(f"✅ Report generated: {output_path}")
        
        # Step 5: Show summary
        print("\n📈 Summary:")
        from collections import Counter
        classifications = Counter(r.classification.value for r in classified_variants)
        
        for classification, count in classifications.most_common():
            print(f"  {classification}: {count}")
        
        print(f"\n🎉 Example completed successfully!")
        print(f"📄 View the report at: {output_path.absolute()}")
        
    except Exception as e:
        print(f"❌ Error during processing: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main() 