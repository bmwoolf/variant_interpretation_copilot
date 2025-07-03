"""
Report generation module for clinical variant interpretation.
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime

from jinja2 import Environment, FileSystemLoader, Template
from rich.console import Console

from .models import ACMGResult, ReportConfig

logger = logging.getLogger(__name__)
console = Console()


class ReportGenerator:
    """Generate clinical variant interpretation reports."""
    
    def __init__(self):
        """Initialize report generator."""
        self.logger = logging.getLogger(__name__)
        
        # Set up Jinja2 environment
        self.env = Environment(
            loader=FileSystemLoader(self._get_template_dir()),
            autoescape=True
        )
    
    def generate(self, results: List[ACMGResult], output_path: Path, config: ReportConfig) -> None:
        """
        Generate report from ACMG results.
        
        Args:
            results: List of ACMG classification results
            output_path: Output file path
            config: Report configuration
        """
        try:
            if config.output_format == "html":
                self._generate_html_report(results, output_path, config)
            elif config.output_format == "json":
                self._generate_json_report(results, output_path, config)
            else:
                raise ValueError(f"Unsupported output format: {config.output_format}")
                
        except Exception as e:
            self.logger.error(f"Report generation failed: {e}")
            raise
    
    def _generate_html_report(self, results: List[ACMGResult], output_path: Path, config: ReportConfig) -> None:
        """Generate HTML report."""
        # Filter results based on configuration
        filtered_results = self._filter_results(results, config)
        
        # Prepare template data
        template_data = {
            'results': filtered_results,
            'config': config,
            'summary': self._generate_summary(filtered_results),
            'generated_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'total_variants': len(results),
            'filtered_variants': len(filtered_results),
        }
        
        # Render template
        template = self.env.get_template('report.html')
        html_content = template.render(**template_data)
        
        # Write to file
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"HTML report generated: {output_path}")
    
    def _generate_json_report(self, results: List[ACMGResult], output_path: Path, config: ReportConfig) -> None:
        """Generate JSON report."""
        # Filter results based on configuration
        filtered_results = self._filter_results(results, config)
        
        # Prepare JSON data
        json_data = {
            'metadata': {
                'generated_at': datetime.now().isoformat(),
                'total_variants': len(results),
                'filtered_variants': len(filtered_results),
                'config': config.dict(),
            },
            'summary': self._generate_summary(filtered_results),
            'variants': [self._result_to_dict(result) for result in filtered_results],
        }
        
        # Write to file
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, default=str)
        
        self.logger.info(f"JSON report generated: {output_path}")
    
    def _filter_results(self, results: List[ACMGResult], config: ReportConfig) -> List[ACMGResult]:
        """Filter results based on configuration."""
        filtered = results
        
        # Filter by quality score
        if config.min_quality is not None:
            filtered = [r for r in filtered if r.variant.qual is None or r.variant.qual >= config.min_quality]
        
        # Filter by gnomAD frequency
        if config.max_gnomad_af is not None:
            filtered = [r for r in filtered if r.variant.gnomad_af is None or r.variant.gnomad_af <= config.max_gnomad_af]
        
        # Filter by phenotype (placeholder for future implementation)
        if config.phenotype_filter:
            # This would typically query HPO database for gene-phenotype associations
            pass
        
        return filtered
    
    def _generate_summary(self, results: List[ACMGResult]) -> Dict[str, Any]:
        """Generate summary statistics."""
        if not results:
            return {
                'total': 0,
                'classifications': {},
                'genes': set(),
                'chromosomes': set(),
            }
        
        # Count classifications
        classifications = {}
        for result in results:
            classification = result.classification.value
            classifications[classification] = classifications.get(classification, 0) + 1
        
        # Collect unique genes and chromosomes
        genes = {result.variant.gene for result in results if result.variant.gene}
        chromosomes = {result.variant.chrom for result in results}
        
        return {
            'total': len(results),
            'classifications': classifications,
            'genes': list(genes),
            'chromosomes': list(chromosomes),
            'gene_count': len(genes),
            'chromosome_count': len(chromosomes),
        }
    
    def _result_to_dict(self, result: ACMGResult) -> Dict[str, Any]:
        """Convert ACMG result to dictionary for JSON serialization."""
        return {
            'variant': {
                'chrom': result.variant.chrom,
                'pos': result.variant.pos,
                'ref': result.variant.ref,
                'alt': result.variant.alt,
                'variant_type': result.variant.variant_type.value,
                'gene': result.variant.gene,
                'transcript': result.variant.transcript,
                'hgvs_c': result.variant.hgvs_c,
                'hgvs_p': result.variant.hgvs_p,
                'impact': result.variant.impact,
                'clinvar_significance': result.variant.clinvar_significance,
                'clinvar_diseases': result.variant.clinvar_diseases,
                'gnomad_af': result.variant.gnomad_af,
                'cadd_score': result.variant.cadd_score,
                'polyphen_score': result.variant.polyphen_score,
                'sift_score': result.variant.sift_score,
            },
            'classification': result.classification.value,
            'criteria': [c.value for c in result.criteria],
            'score': result.score,
            'reasoning': result.reasoning,
            'confidence': result.confidence,
        }
    
    def _get_template_dir(self) -> Path:
        """Get template directory path."""
        # Look for templates in the package directory
        package_dir = Path(__file__).parent
        template_dir = package_dir / 'templates'
        
        if not template_dir.exists():
            # Create default template
            self._create_default_template(template_dir)
        
        return template_dir
    
    def _create_default_template(self, template_dir: Path) -> None:
        """Create default HTML template."""
        template_dir.mkdir(exist_ok=True)
        
        html_template = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Variant Interpretation Report</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }
        .header h1 {
            margin: 0;
            font-size: 2.5em;
            font-weight: 300;
        }
        .header p {
            margin: 10px 0 0 0;
            opacity: 0.9;
        }
        .summary {
            padding: 20px;
            background: #f8f9fa;
            border-bottom: 1px solid #dee2e6;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }
        .summary-card {
            background: white;
            padding: 15px;
            border-radius: 6px;
            text-align: center;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        .summary-card h3 {
            margin: 0 0 10px 0;
            color: #495057;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .summary-card .value {
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }
        .variants {
            padding: 20px;
        }
        .variant-card {
            background: white;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            margin-bottom: 20px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }
        .variant-header {
            padding: 15px 20px;
            border-bottom: 1px solid #dee2e6;
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: #f8f9fa;
        }
        .variant-info {
            flex: 1;
        }
        .variant-location {
            font-weight: bold;
            color: #495057;
            font-size: 1.1em;
        }
        .variant-gene {
            color: #6c757d;
            font-size: 0.9em;
            margin-top: 2px;
        }
        .classification-badge {
            padding: 6px 12px;
            border-radius: 20px;
            font-size: 0.8em;
            font-weight: bold;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .classification-pathogenic {
            background: #dc3545;
            color: white;
        }
        .classification-likely-pathogenic {
            background: #fd7e14;
            color: white;
        }
        .classification-uncertain {
            background: #ffc107;
            color: #212529;
        }
        .classification-likely-benign {
            background: #20c997;
            color: white;
        }
        .classification-benign {
            background: #28a745;
            color: white;
        }
        .variant-details {
            padding: 20px;
        }
        .detail-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin-bottom: 15px;
        }
        .detail-item {
            display: flex;
            justify-content: space-between;
            padding: 8px 0;
            border-bottom: 1px solid #f1f3f4;
        }
        .detail-label {
            font-weight: 500;
            color: #495057;
        }
        .detail-value {
            color: #6c757d;
        }
        .criteria-section {
            margin-top: 15px;
            padding-top: 15px;
            border-top: 2px solid #e9ecef;
        }
        .criteria-title {
            font-weight: bold;
            color: #495057;
            margin-bottom: 10px;
        }
        .criteria-list {
            list-style: none;
            padding: 0;
            margin: 0;
        }
        .criteria-item {
            background: #f8f9fa;
            padding: 8px 12px;
            margin-bottom: 5px;
            border-radius: 4px;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
        }
        .reasoning {
            margin-top: 15px;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 6px;
            font-style: italic;
            color: #495057;
        }
        .footer {
            padding: 20px;
            text-align: center;
            color: #6c757d;
            font-size: 0.9em;
            border-top: 1px solid #dee2e6;
        }
        @media (max-width: 768px) {
            .variant-header {
                flex-direction: column;
                align-items: flex-start;
                gap: 10px;
            }
            .summary-grid {
                grid-template-columns: 1fr;
            }
            .detail-grid {
                grid-template-columns: 1fr;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Clinical Variant Interpretation Report</h1>
            <p>Generated on {{ generated_at }}</p>
        </div>
        
        <div class="summary">
            <h2>Summary</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>Total Variants</h3>
                    <div class="value">{{ summary.total }}</div>
                </div>
                <div class="summary-card">
                    <h3>Genes</h3>
                    <div class="value">{{ summary.gene_count }}</div>
                </div>
                <div class="summary-card">
                    <h3>Chromosomes</h3>
                    <div class="value">{{ summary.chromosome_count }}</div>
                </div>
                <div class="summary-card">
                    <h3>Filtered</h3>
                    <div class="value">{{ filtered_variants }}/{{ total_variants }}</div>
                </div>
            </div>
            
            {% if summary.classifications %}
            <h3>Classifications</h3>
            <div class="summary-grid">
                {% for classification, count in summary.classifications.items() %}
                <div class="summary-card">
                    <h3>{{ classification }}</h3>
                    <div class="value">{{ count }}</div>
                </div>
                {% endfor %}
            </div>
            {% endif %}
        </div>
        
        <div class="variants">
            <h2>Variants</h2>
            {% for result in results %}
            <div class="variant-card">
                <div class="variant-header">
                    <div class="variant-info">
                        <div class="variant-location">{{ result.variant.chrom }}:{{ result.variant.pos }} {{ result.variant.ref }}>{{ result.variant.alt }}</div>
                        {% if result.variant.gene %}
                        <div class="variant-gene">{{ result.variant.gene }}</div>
                        {% endif %}
                    </div>
                    <div class="classification-badge classification-{{ result.classification.value.lower().replace(' ', '-') }}">
                        {{ result.classification.value }}
                    </div>
                </div>
                
                <div class="variant-details">
                    <div class="detail-grid">
                        <div class="detail-item">
                            <span class="detail-label">Type:</span>
                            <span class="detail-value">{{ result.variant.variant_type.value }}</span>
                        </div>
                        {% if result.variant.hgvs_c %}
                        <div class="detail-item">
                            <span class="detail-label">cDNA:</span>
                            <span class="detail-value">{{ result.variant.hgvs_c }}</span>
                        </div>
                        {% endif %}
                        {% if result.variant.hgvs_p %}
                        <div class="detail-item">
                            <span class="detail-label">Protein:</span>
                            <span class="detail-value">{{ result.variant.hgvs_p }}</span>
                        </div>
                        {% endif %}
                        {% if result.variant.impact %}
                        <div class="detail-item">
                            <span class="detail-label">Impact:</span>
                            <span class="detail-value">{{ result.variant.impact }}</span>
                        </div>
                        {% endif %}
                        {% if result.variant.clinvar_significance %}
                        <div class="detail-item">
                            <span class="detail-label">ClinVar:</span>
                            <span class="detail-value">{{ result.variant.clinvar_significance }}</span>
                        </div>
                        {% endif %}
                        {% if result.variant.gnomad_af is not none %}
                        <div class="detail-item">
                            <span class="detail-label">gnomAD AF:</span>
                            <span class="detail-value">{{ "%.6f"|format(result.variant.gnomad_af) }}</span>
                        </div>
                        {% endif %}
                        {% if result.variant.cadd_score is not none %}
                        <div class="detail-item">
                            <span class="detail-label">CADD:</span>
                            <span class="detail-value">{{ "%.2f"|format(result.variant.cadd_score) }}</span>
                        </div>
                        {% endif %}
                        <div class="detail-item">
                            <span class="detail-label">ACMG Score:</span>
                            <span class="detail-value">{{ result.score }}</span>
                        </div>
                        <div class="detail-item">
                            <span class="detail-label">Confidence:</span>
                            <span class="detail-value">{{ "%.1f"|format(result.confidence * 100) }}%</span>
                        </div>
                    </div>
                    
                    {% if result.criteria %}
                    <div class="criteria-section">
                        <div class="criteria-title">Applied ACMG Criteria:</div>
                        <ul class="criteria-list">
                            {% for criterion in result.criteria %}
                            <li class="criteria-item">{{ criterion.value }}</li>
                            {% endfor %}
                        </ul>
                    </div>
                    {% endif %}
                    
                    {% if result.reasoning %}
                    <div class="reasoning">
                        {{ result.reasoning }}
                    </div>
                    {% endif %}
                </div>
            </div>
            {% endfor %}
        </div>
        
        <div class="footer">
            <p>Generated by VCF Copilot - Clinical Genomics Variant Interpretation Tool</p>
        </div>
    </div>
</body>
</html>'''
        
        template_path = template_dir / 'report.html'
        with open(template_path, 'w', encoding='utf-8') as f:
            f.write(html_template)
        
        self.logger.info(f"Default HTML template created: {template_path}") 