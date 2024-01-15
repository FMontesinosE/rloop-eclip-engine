# Development Guide

This document provides instructions for developers working on rloop-eclip-engine.

## Setup Development Environment

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone repository
git clone https://github.com/yourusername/rloop-eclip-engine.git
cd rloop-eclip-engine

# Create virtual environment and install with dev dependencies
uv venv
source .venv/bin/activate
uv sync --all-extras
```

## Code Quality Standards

### Linting and Formatting

```bash
# Format code (applies changes)
ruff format src/ tests/

# Check for lint errors
ruff check src/ tests/

# Auto-fix lint issues where possible
ruff check src/ tests/ --fix
```

### Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov=src --cov-report=html

# Run only fast tests (exclude slow integration tests)
pytest tests/ -v -m "not slow"

# Run benchmarks only
pytest tests/ --benchmark-only
```

## Project Structure

```
rloop-eclip-engine/
├── src/rloop_eclip/
│   ├── __init__.py          # Package initialization
│   ├── cli.py               # Click CLI interface
│   ├── core.py              # Chromosome normalization utilities
│   ├── engine.py            # Main enrichment engine
│   └── schemas.py           # Pydantic validation models
├── tests/
│   ├── __init__.py
│   └── test_engine.py       # Unit and integration tests
├── config/
│   └── example_config.json  # Sample configuration
├── scripts/
│   └── quickstart.sh        # Quick start demo script
├── pyproject.toml           # Package configuration
├── Dockerfile               # Container build definition
└── README.md                # User documentation
```

## Architecture Overview

### Data Flow

1. **Input Loading** (`core.py`):
   - Parse BED files into Polars LazyFrames
   - Normalize chromosome names (UCSC ↔ Ensembl)
   - Filter excluded regions and invalid intervals

2. **Overlap Calculation** (`engine.py`):
   - Join intervals on matching chromosomes
   - Calculate base-pair overlaps using vectorized operations
   - Return observed overlap statistic

3. **Permutation Testing** (`engine.py`):
   - Stratified shuffling by chromosome
   - Preserve interval widths during randomization
   - Generate null distribution via Monte Carlo sampling
   - Compute empirical p-value with epsilon smoothing

4. **Result Output** (`cli.py`):
   - Serialize EnrichmentResult to Parquet or JSON
   - Display formatted table in terminal via Rich

### Key Design Decisions

- **Polars over Pandas**: Lazy evaluation enables out-of-core processing of large genomic datasets without OOM errors
- **PyTorch for Permutations**: Tensor operations reduce GC overhead during 10k+ iteration resampling
- **Pydantic v2 Validation**: Strict schema enforcement catches configuration errors before expensive computation
- **Loguru Logging**: Structured logs with timestamps for production debugging

## Performance Optimization

### Memory Profiling

```bash
# Install memory_profiler
uv pip install memory_profiler

# Profile specific function
python -m memory_profiler scripts/profile_overlap.py
```

### Benchmark Comparison

```bash
# Run benchmarks and save baseline
pytest tests/test_engine.py --benchmark-only --benchmark-save=baseline

# After optimization, compare
pytest tests/test_engine.py --benchmark-only --benchmark-compare=baseline
```

## Contributing Guidelines

1. **Branch Naming**: Use conventional prefixes (`feat/`, `fix/`, `refactor/`)
2. **Commit Messages**: Follow Conventional Commits specification
3. **Pre-commit Checks**: Run `ruff format` and `pytest` before committing
4. **Code Reviews**: All changes require review before merge to main

## Debugging Tips

### Enable Verbose Logging

```bash
rloop-eclip --verbose \
    --eclip-bed data/test.bed \
    --drip-bed data/rloops.bed \
    --genome-sizes hg38.chrom.sizes \
    --output debug_results.parquet
```

### Interactive Debugging

```python
# Add breakpoint in code
import pdb; pdb.set_trace()

# Or use ipdb for better interface
import ipdb; ipdb.set_trace()
```

### Common Issues

**Issue**: "No valid chromosomes after filtering"
- **Cause**: Chromosome naming mismatch between BED files
- **Fix**: Check if one file uses UCSC (chr1) and other uses Ensembl (1)

**Issue**: "Memory allocation failed"
- **Cause**: Dataset too large for available RAM
- **Fix**: Reduce `--iterations` or increase swap space

**Issue**: "Empirical p-value is exactly 0.0001"
- **Cause**: Not enough permutations to resolve small p-values
- **Fix**: Increase `--iterations` to 100000 for publication-quality results
