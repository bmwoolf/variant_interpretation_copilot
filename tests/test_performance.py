"""
Performance benchmarking tests for VCF Copilot.

This module provides comprehensive benchmarking to measure:
- VCF parsing speed and throughput
- Memory usage during processing
- Annotation performance
- ACMG classification speed
- End-to-end pipeline performance

Results can be compared against future Rust implementation.
"""

import time
import psutil
import pytest
import tempfile
import subprocess
from pathlib import Path
from typing import List, Dict, Any
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
import statistics

from vcf_copilot.parser import VCFParser
from vcf_copilot.annotators.engine import AnnotationEngine
from vcf_copilot.scoring import ACMGClassifier
from vcf_copilot.models import Variant, VariantType


@dataclass
class BenchmarkResult:
    """Result of a performance benchmark."""
    operation: str
    duration: float
    memory_peak: float
    memory_avg: float
    throughput: float  # operations per second
    variant_count: int
    metadata: Dict[str, Any]


class PerformanceBenchmark:
    """Performance benchmarking suite for VCF Copilot."""
    
    def __init__(self):
        self.parser = VCFParser()
        self.annotator = AnnotationEngine()
        self.classifier = ACMGClassifier()
        self.process = psutil.Process()
    
    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB."""
        return self.process.memory_info().rss / 1024 / 1024
    
    def _measure_memory(self, func, *args, **kwargs) -> tuple:
        """Measure memory usage during function execution."""
        memory_samples = []
        
        def memory_monitor():
            while True:
                memory_samples.append(self._get_memory_usage())
                time.sleep(0.01)  # Sample every 10ms
        
        import threading
        monitor_thread = threading.Thread(target=memory_monitor, daemon=True)
        monitor_thread.start()
        
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        
        # Stop monitoring
        monitor_thread.join(timeout=0.1)
        
        duration = end_time - start_time
        memory_peak = max(memory_samples) if memory_samples else 0
        memory_avg = statistics.mean(memory_samples) if memory_samples else 0
        
        return result, duration, memory_peak, memory_avg
    
    def benchmark_vcf_parsing(self, vcf_path: Path, iterations: int = 3) -> BenchmarkResult:
        """Benchmark VCF parsing performance."""
        print(f"Benchmarking VCF parsing: {vcf_path}")
        
        durations = []
        memory_peaks = []
        memory_avgs = []
        variant_counts = []
        
        for i in range(iterations):
            print(f"  Iteration {i+1}/{iterations}")
            
            # Clear memory before each iteration
            import gc
            gc.collect()
            
            result, duration, memory_peak, memory_avg = self._measure_memory(
                self.parser.parse, vcf_path
            )
            
            durations.append(duration)
            memory_peaks.append(memory_peak)
            memory_avgs.append(memory_avg)
            variant_counts.append(len(result))
        
        avg_duration = statistics.mean(durations)
        avg_memory_peak = statistics.mean(memory_peaks)
        avg_memory_avg = statistics.mean(memory_avgs)
        avg_variant_count = statistics.mean(variant_counts)
        
        throughput = avg_variant_count / avg_duration if avg_duration > 0 else 0
        
        return BenchmarkResult(
            operation="VCF Parsing",
            duration=avg_duration,
            memory_peak=avg_memory_peak,
            memory_avg=avg_memory_avg,
            throughput=throughput,
            variant_count=int(avg_variant_count),
            metadata={
                "iterations": iterations,
                "file_size_mb": vcf_path.stat().st_size / 1024 / 1024,
                "duration_std": statistics.stdev(durations),
                "memory_peak_std": statistics.stdev(memory_peaks),
            }
        )
    
    def benchmark_annotation(self, variants: List[Variant], iterations: int = 3) -> BenchmarkResult:
        """Benchmark annotation performance."""
        print(f"Benchmarking annotation for {len(variants)} variants")
        
        durations = []
        memory_peaks = []
        memory_avgs = []
        
        for i in range(iterations):
            print(f"  Iteration {i+1}/{iterations}")
            
            # Clear memory before each iteration
            import gc
            gc.collect()
            
            result, duration, memory_peak, memory_avg = self._measure_memory(
                self.annotator.annotate_batch, variants
            )
            
            durations.append(duration)
            memory_peaks.append(memory_peak)
            memory_avgs.append(memory_avg)
        
        avg_duration = statistics.mean(durations)
        avg_memory_peak = statistics.mean(memory_peaks)
        avg_memory_avg = statistics.mean(memory_avgs)
        
        throughput = len(variants) / avg_duration if avg_duration > 0 else 0
        
        return BenchmarkResult(
            operation="Annotation",
            duration=avg_duration,
            memory_peak=avg_memory_peak,
            memory_avg=avg_memory_avg,
            throughput=throughput,
            variant_count=len(variants),
            metadata={
                "iterations": iterations,
                "duration_std": statistics.stdev(durations),
                "memory_peak_std": statistics.stdev(memory_peaks),
            }
        )
    
    def benchmark_classification(self, variants: List[Variant], iterations: int = 3) -> BenchmarkResult:
        """Benchmark ACMG classification performance."""
        print(f"Benchmarking classification for {len(variants)} variants")
        
        durations = []
        memory_peaks = []
        memory_avgs = []
        
        for i in range(iterations):
            print(f"  Iteration {i+1}/{iterations}")
            
            # Clear memory before each iteration
            import gc
            gc.collect()
            
            def classify_variants():
                return [self.classifier.classify(variant) for variant in variants]
            
            result, duration, memory_peak, memory_avg = self._measure_memory(classify_variants)
            
            durations.append(duration)
            memory_peaks.append(memory_peak)
            memory_avgs.append(memory_avg)
        
        avg_duration = statistics.mean(durations)
        avg_memory_peak = statistics.mean(memory_peaks)
        avg_memory_avg = statistics.mean(memory_avgs)
        
        throughput = len(variants) / avg_duration if avg_duration > 0 else 0
        
        return BenchmarkResult(
            operation="ACMG Classification",
            duration=avg_duration,
            memory_peak=avg_memory_peak,
            memory_avg=avg_memory_avg,
            throughput=throughput,
            variant_count=len(variants),
            metadata={
                "iterations": iterations,
                "duration_std": statistics.stdev(durations),
                "memory_peak_std": statistics.stdev(memory_peaks),
            }
        )
    
    def benchmark_end_to_end(self, vcf_path: Path, iterations: int = 3) -> BenchmarkResult:
        """Benchmark complete end-to-end pipeline."""
        print(f"Benchmarking end-to-end pipeline: {vcf_path}")
        
        durations = []
        memory_peaks = []
        memory_avgs = []
        variant_counts = []
        
        for i in range(iterations):
            print(f"  Iteration {i+1}/{iterations}")
            
            # Clear memory before each iteration
            import gc
            gc.collect()
            
            def run_pipeline():
                # Parse
                variants = self.parser.parse(vcf_path)
                # Annotate
                annotated_variants = self.annotator.annotate_batch(variants)
                # Classify
                classified_variants = [self.classifier.classify(v) for v in annotated_variants]
                return classified_variants
            
            result, duration, memory_peak, memory_avg = self._measure_memory(run_pipeline)
            
            durations.append(duration)
            memory_peaks.append(memory_peak)
            memory_avgs.append(memory_avg)
            variant_counts.append(len(result))
        
        avg_duration = statistics.mean(durations)
        avg_memory_peak = statistics.mean(memory_peaks)
        avg_memory_avg = statistics.mean(memory_avgs)
        avg_variant_count = statistics.mean(variant_counts)
        
        throughput = avg_variant_count / avg_duration if avg_duration > 0 else 0
        
        return BenchmarkResult(
            operation="End-to-End Pipeline",
            duration=avg_duration,
            memory_peak=avg_memory_peak,
            memory_avg=avg_memory_avg,
            throughput=throughput,
            variant_count=int(avg_variant_count),
            metadata={
                "iterations": iterations,
                "file_size_mb": vcf_path.stat().st_size / 1024 / 1024,
                "duration_std": statistics.stdev(durations),
                "memory_peak_std": statistics.stdev(memory_peaks),
            }
        )
    
    def generate_large_vcf(self, variant_count: int = 10000, output_path: Path = None) -> Path:
        """Generate a large VCF file for benchmarking."""
        if output_path is None:
            output_path = Path(f"data/benchmark_{variant_count}_variants.vcf")
        
        output_path.parent.mkdir(exist_ok=True)
        
        with open(output_path, 'w') as f:
            # Write header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
            f.write("##INFO=<ID=Gene_Name,Number=1,Type=String,Description=\"Gene name\">\n")
            f.write("##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Variant impact\">\n")
            f.write("##INFO=<ID=CADD_RAW,Number=1,Type=Float,Description=\"CADD raw score\">\n")
            f.write("##INFO=<ID=PolyPhen,Number=1,Type=String,Description=\"PolyPhen prediction\">\n")
            f.write("##INFO=<ID=SIFT,Number=1,Type=String,Description=\"SIFT prediction\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n")
            
            # Generate variants
            import random
            
            chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
            genes = ["BRCA1", "BRCA2", "TP53", "PTEN", "APC", "RB1", "VHL", "ATM", "CHEK2", "PALB2"]
            impacts = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
            
            for i in range(variant_count):
                chrom = random.choice(chromosomes)
                pos = random.randint(1, 100000000)
                ref = random.choice(["A", "T", "G", "C"])
                alt = random.choice(["A", "T", "G", "C"])
                gene = random.choice(genes)
                impact = random.choice(impacts)
                cadd_score = random.uniform(0, 30)
                
                info = f"Gene_Name={gene};IMPACT={impact};CADD_RAW={cadd_score:.2f};PolyPhen=probably_damaging:0.95;SIFT=deleterious:0.01"
                
                f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\t{info}\tGT:AD:DP:GQ:PL\t0/1:45,55:100:99:1000,0,2000\n")
        
        return output_path
    
    def run_comprehensive_benchmark(self, vcf_path: Path = None) -> List[BenchmarkResult]:
        """Run comprehensive benchmark suite."""
        if vcf_path is None:
            vcf_path = Path("data/test.vcf")
        
        if not vcf_path.exists():
            pytest.skip(f"VCF file not found: {vcf_path}")
        
        results = []
        
        # Benchmark parsing
        results.append(self.benchmark_vcf_parsing(vcf_path))
        
        # Parse variants for other benchmarks
        variants = self.parser.parse(vcf_path)
        
        # Benchmark annotation
        results.append(self.benchmark_annotation(variants))
        
        # Benchmark classification
        results.append(self.benchmark_classification(variants))
        
        # Benchmark end-to-end
        results.append(self.benchmark_end_to_end(vcf_path))
        
        return results
    
    def print_benchmark_results(self, results: List[BenchmarkResult]):
        """Print benchmark results in a formatted table."""
        print("\n" + "="*80)
        print("PERFORMANCE BENCHMARK RESULTS")
        print("="*80)
        
        print(f"{'Operation':<25} {'Duration (s)':<15} {'Throughput':<15} {'Memory Peak (MB)':<18} {'Variants':<10}")
        print("-"*80)
        
        for result in results:
            print(f"{result.operation:<25} {result.duration:<15.3f} {result.throughput:<15.1f} {result.memory_peak:<18.1f} {result.variant_count:<10}")
        
        print("\n" + "="*80)
        print("DETAILED RESULTS")
        print("="*80)
        
        for result in results:
            print(f"\n{result.operation}:")
            print(f"  Duration: {result.duration:.3f}s ± {result.metadata.get('duration_std', 0):.3f}s")
            print(f"  Throughput: {result.throughput:.1f} variants/second")
            print(f"  Memory Peak: {result.memory_peak:.1f} MB ± {result.metadata.get('memory_peak_std', 0):.1f} MB")
            print(f"  Memory Average: {result.memory_avg:.1f} MB")
            print(f"  Variant Count: {result.variant_count}")
            
            if 'file_size_mb' in result.metadata:
                print(f"  File Size: {result.metadata['file_size_mb']:.2f} MB")
                mb_per_second = result.metadata['file_size_mb'] / result.duration
                print(f"  Processing Speed: {mb_per_second:.2f} MB/s")


class TestPerformance:
    """Performance benchmark tests."""
    
    def setup_method(self):
        """Set up benchmark suite."""
        self.benchmark = PerformanceBenchmark()
    
    def test_benchmark_small_vcf(self):
        """Benchmark performance on small VCF file."""
        vcf_path = Path("data/test.vcf")
        if not vcf_path.exists():
            pytest.skip("Test VCF file not found")
        
        results = self.benchmark.run_comprehensive_benchmark(vcf_path)
        self.benchmark.print_benchmark_results(results)
        
        # Basic assertions to ensure reasonable performance
        for result in results:
            assert result.duration > 0, f"Duration should be positive for {result.operation}"
            assert result.throughput > 0, f"Throughput should be positive for {result.operation}"
            assert result.memory_peak > 0, f"Memory usage should be positive for {result.operation}"
    
    def test_benchmark_large_vcf(self):
        """Benchmark performance on large VCF file."""
        # Generate large VCF if it doesn't exist
        large_vcf_path = Path("data/benchmark_10000_variants.vcf")
        if not large_vcf_path.exists():
            print("Generating large VCF file for benchmarking...")
            large_vcf_path = self.benchmark.generate_large_vcf(10000, large_vcf_path)
        
        results = self.benchmark.run_comprehensive_benchmark(large_vcf_path)
        self.benchmark.print_benchmark_results(results)
        
        # Performance assertions for large file
        parsing_result = next(r for r in results if r.operation == "VCF Parsing")
        assert parsing_result.throughput > 1000, f"Should parse at least 1000 variants/second, got {parsing_result.throughput}"
        
        e2e_result = next(r for r in results if r.operation == "End-to-End Pipeline")
        assert e2e_result.throughput > 100, f"Should process at least 100 variants/second end-to-end, got {e2e_result.throughput}"
    
    def test_memory_efficiency(self):
        """Test memory efficiency during processing."""
        vcf_path = Path("data/test.vcf")
        if not vcf_path.exists():
            pytest.skip("Test VCF file not found")
        
        # Measure memory before
        initial_memory = self.benchmark._get_memory_usage()
        
        # Run pipeline
        variants = self.benchmark.parser.parse(vcf_path)
        annotated_variants = self.benchmark.annotator.annotate_batch(variants)
        classified_variants = [self.benchmark.classifier.classify(v) for v in annotated_variants]
        
        # Measure memory after
        final_memory = self.benchmark._get_memory_usage()
        memory_increase = final_memory - initial_memory
        
        print(f"Memory usage: {initial_memory:.1f} MB -> {final_memory:.1f} MB (+{memory_increase:.1f} MB)")
        
        # Memory should not increase excessively
        assert memory_increase < 100, f"Memory increase should be less than 100 MB, got {memory_increase:.1f} MB"
    
    def test_scalability(self):
        """Test performance scalability with different file sizes."""
        sizes = [100, 1000, 5000]
        results = {}
        
        for size in sizes:
            vcf_path = Path(f"data/benchmark_{size}_variants.vcf")
            if not vcf_path.exists():
                vcf_path = self.benchmark.generate_large_vcf(size, vcf_path)
            
            # Benchmark parsing only for scalability test
            result = self.benchmark.benchmark_vcf_parsing(vcf_path, iterations=2)
            results[size] = result
        
        # Print scalability results
        print("\n" + "="*60)
        print("SCALABILITY TEST RESULTS")
        print("="*60)
        print(f"{'Variants':<10} {'Duration (s)':<15} {'Throughput':<15} {'Memory (MB)':<15}")
        print("-"*60)
        
        for size, result in results.items():
            print(f"{size:<10} {result.duration:<15.3f} {result.throughput:<15.1f} {result.memory_peak:<15.1f}")
        
        # Check that throughput doesn't degrade too much with size
        small_throughput = results[100].throughput
        large_throughput = results[5000].throughput
        
        # Throughput should not degrade more than 50%
        degradation = (small_throughput - large_throughput) / small_throughput
        assert degradation < 0.5, f"Throughput degradation should be less than 50%, got {degradation:.1%}"


def run_benchmark_cli():
    """CLI function to run benchmarks."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Run VCF Copilot performance benchmarks")
    parser.add_argument("--vcf", type=Path, default="data/test.vcf", help="VCF file to benchmark")
    parser.add_argument("--large", action="store_true", help="Generate and test large VCF file")
    parser.add_argument("--scalability", action="store_true", help="Run scalability tests")
    
    args = parser.parse_args()
    
    benchmark = PerformanceBenchmark()
    
    if args.large:
        print("Generating large VCF file...")
        large_vcf = benchmark.generate_large_vcf(10000)
        results = benchmark.run_comprehensive_benchmark(large_vcf)
    elif args.scalability:
        print("Running scalability tests...")
        test = TestPerformance()
        test.setup_method()
        test.test_scalability()
        return
    else:
        results = benchmark.run_comprehensive_benchmark(args.vcf)
    
    benchmark.print_benchmark_results(results)


if __name__ == "__main__":
    run_benchmark_cli() 