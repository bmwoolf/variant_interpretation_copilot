.PHONY: help install install-dev test lint format clean example run benchmark

help: ## Show this help message
	@echo "VCF Copilot - Clinical Genomics Variant Interpretation Tool"
	@echo ""
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install: ## Install the package in development mode
	pip install -e .

install-dev: ## Install the package with development dependencies
	pip install -e ".[dev]"

test: ## Run tests
	pytest tests/ -v --cov=vcf_copilot --cov-report=term-missing

lint: ## Run linting checks
	flake8 vcf_copilot/ tests/
	isort --check-only vcf_copilot/ tests/
	black --check vcf_copilot/ tests/

format: ## Format code
	isort vcf_copilot/ tests/
	black vcf_copilot/ tests/

clean: ## Clean up generated files
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} +
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -f example_report.html
	rm -f report.html
	rm -f variants.json

example: ## Run the example script
	python example.py

run: ## Run VCF Copilot on test data
	python -m vcf_copilot data/test.vcf --output report.html --verbose

validate: ## Validate the test VCF file
	python -m vcf_copilot validate data/test.vcf

version: ## Show version information
	python -m vcf_copilot version

setup: install-dev ## Set up development environment
	@echo "Development environment set up successfully!"
	@echo "Run 'make example' to test the installation."

benchmark: ## Run performance benchmarks
	pytest -s tests/test_performance.py 