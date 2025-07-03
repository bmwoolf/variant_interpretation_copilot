"""
Main CLI application for VCF Copilot.
"""

import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

from .models import ReportConfig
from .parser import VCFParser
from .annotators.engine import AnnotationEngine
from .scoring import ACMGClassifier
from .report import ReportGenerator

app = typer.Typer(
    name="vcf-copilot",
    help="Clinical Genomics Copilot for VCF variant interpretation",
    add_completion=False,
)
console = Console()


@app.command()
def main(
    input_file: Path = typer.Argument(
        ...,
        help="Input VCF file",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        "report.html",
        "--output",
        "-o",
        help="Output file path",
    ),
    format: str = typer.Option(
        "html",
        "--format",
        "-f",
        help="Output format (html/json)",
        case_sensitive=False,
    ),
    phenotype: Optional[str] = typer.Option(
        None,
        "--phenotype",
        "-p",
        help="HPO phenotype filter (e.g., HP:0001250)",
    ),
    min_quality: Optional[float] = typer.Option(
        None,
        "--min-quality",
        help="Minimum quality score filter",
    ),
    max_gnomad_af: Optional[float] = typer.Option(
        0.01,
        "--max-gnomad-af",
        help="Maximum gnomAD allele frequency for reporting",
    ),
    include_evidence: bool = typer.Option(
        True,
        "--include-evidence/--no-evidence",
        help="Include evidence details in report",
    ),
    include_recommendations: bool = typer.Option(
        True,
        "--include-recommendations/--no-recommendations",
        help="Include clinical recommendations in report",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Verbose output",
    ),
) -> None:
    """
    Process VCF file and generate clinical variant interpretation report.
    
    Example:
        vcf-copilot input.vcf --output report.html --phenotype HP:0001250
    """
    try:
        # Validate inputs
        if format.lower() not in ["html", "json"]:
            console.print(f"[red]Error: Unsupported format '{format}'. Use 'html' or 'json'.[/red]")
            sys.exit(1)
        
        # Create report configuration
        config = ReportConfig(
            output_format=format.lower(),
            include_evidence=include_evidence,
            include_recommendations=include_recommendations,
            phenotype_filter=phenotype,
            min_quality=min_quality,
            max_gnomad_af=max_gnomad_af,
        )
        
        console.print(f"[bold blue]🧬 Clinical Genomics Copilot[/bold blue]")
        console.print(f"Processing: [bold]{input_file}[/bold]")
        console.print(f"Output: [bold]{output}[/bold]")
        console.print()
        
        # Process VCF file
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
            disable=verbose,
        ) as progress:
            # Parse VCF
            task = progress.add_task("Parsing VCF file...", total=None)
            parser = VCFParser()
            variants = parser.parse(input_file)
            progress.update(task, description=f"Parsed {len(variants)} variants")
            
            if verbose:
                console.print(f"Found {len(variants)} variants")
            
            # Annotate variants
            task = progress.add_task("Annotating variants...", total=len(variants))
            annotator = AnnotationEngine()
            annotated_variants = []
            
            for variant in variants:
                annotated_variant = annotator.annotate(variant)
                annotated_variants.append(annotated_variant)
                progress.advance(task)
            
            # Classify variants
            task = progress.add_task("Classifying variants...", total=len(annotated_variants))
            classifier = ACMGClassifier()
            classified_variants = []
            
            for variant in annotated_variants:
                classified_variant = classifier.classify(variant)
                classified_variants.append(classified_variant)
                progress.advance(task)
            
            # Generate report
            task = progress.add_task("Generating report...", total=None)
            report_generator = ReportGenerator()
            report_generator.generate(classified_variants, output, config)
            progress.update(task, description="Report generated successfully")
        
        # Summary
        console.print()
        console.print("[bold green]✅ Processing complete![/bold green]")
        
        # Show summary table
        summary_table = Table(title="Variant Summary")
        summary_table.add_column("Classification", style="cyan")
        summary_table.add_column("Count", justify="right", style="green")
        
        from collections import Counter
        classifications = Counter(v.variant.acmg_classification for v in classified_variants)
        
        for classification, count in classifications.most_common():
            summary_table.add_row(classification or "Unclassified", str(count))
        
        console.print(summary_table)
        
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        if verbose:
            import traceback
            console.print(traceback.format_exc())
        sys.exit(1)


@app.command()
def version() -> None:
    """Show version information."""
    from . import __version__
    console.print(f"VCF Copilot version {__version__}")


@app.command()
def validate(
    input_file: Path = typer.Argument(
        ...,
        help="Input VCF file to validate",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
) -> None:
    """Validate VCF file format and content."""
    try:
        console.print(f"Validating VCF file: [bold]{input_file}[/bold]")
        
        parser = VCFParser()
        validation_result = parser.validate(input_file)
        
        if validation_result.is_valid:
            console.print("[green]✅ VCF file is valid[/green]")
            console.print(f"Variants found: {validation_result.variant_count}")
            console.print(f"Sample count: {validation_result.sample_count}")
        else:
            console.print("[red]❌ VCF file has issues:[/red]")
            for error in validation_result.errors:
                console.print(f"  - {error}")
            sys.exit(1)
            
    except Exception as e:
        console.print(f"[red]Error during validation: {e}[/red]")
        sys.exit(1)


if __name__ == "__main__":
    app() 