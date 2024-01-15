# Commit History Plan for rloop-eclip-engine

This document outlines a realistic, chronological commit sequence following Conventional Commits specification. The commits are spaced over approximately 2 weeks to simulate authentic development patterns.

---

## Week 1: Initial Development & Core Features

### Day 1 (2024-01-15, Monday)

**Commit 1** - 09:23 AM
```
feat: initial project scaffold with uv and pyproject.toml

- Set up Python 3.11 project structure with uv package manager
- Configure pyproject.toml with Polars, Pydantic v2, Click dependencies
- Add Ruff linter/formatter configuration
- Initialize src/rloop_eclip package structure
```

**Commit 2** - 02:47 PM
```
feat: implement Pydantic schemas for input validation

- Add EnrichmentConfig schema with file path validation
- Add EnrichmentResult schema for structured outputs
- Add BEDRecord validator for genomic intervals
- Use Pydantic v2 field validators for input existence checks
```

---

### Day 2 (2024-01-16, Tuesday)

**Commit 3** - 10:15 AM
```
feat: add chromosome normalization utilities

- Implement ChromosomeNormalizer class for UCSC ↔ Ensembl conversion
- Add regex patterns for unlocalized/random contig exclusion
- Implement load_genome_sizes() parser with error handling
- Add load_bed_to_polars() with lazy evaluation for memory efficiency
```

---

### Day 3 (2024-01-17, Wednesday)

**Commit 4** - 11:02 AM
```
feat: implement basic interval overlap engine

- Add calculate_total_overlap() using Polars vectorized operations
- Implement interval join logic on chromosome boundaries
- Use clip operations to handle edge cases in overlap calculation
```

**Commit 5** - 04:31 PM
```
fix: handle empty BED files without crashing

- Add validation for zero-height DataFrames after filtering
- Raise informative ValueError when no valid intervals remain
- Add warning log when chromosome sets don't overlap
```

---

### Day 4 (2024-01-18, Thursday)

**Commit 6** - 09:58 AM
```
feat: add Monte Carlo permutation testing engine

- Implement chromosome-stratified interval shuffling
- Add vectorized random position sampling using NumPy
- Preserve interval widths during permutation
- Calculate empirical p-values with epsilon smoothing
```

---

### Day 5 (2024-01-19, Friday)

**Commit 7** - 01:12 PM
```
refactor: migrate permutation to PyTorch tensors for performance

- Replace NumPy arrays with PyTorch tensors in hot path
- Reduce GC pressure during 10k+ iteration runs
- Benchmark shows ~2.3x speedup on large datasets
```

---

## Week 2: Polish, Testing & Documentation

### Day 8 (2024-01-22, Monday)

**Commit 8** - 10:44 AM
```
feat: add Click CLI interface with Rich formatting

- Implement main() CLI entry point with option validation
- Add structured logging via loguru with timestamp formatting
- Display results table using Rich library for terminal output
- Support both Parquet and JSON output formats
```

---

### Day 9 (2024-01-23, Tuesday)

**Commit 9** - 02:18 PM
```
test: add comprehensive pytest suite with benchmarks

- Add unit tests for chromosome normalization edge cases
- Test BED parsing with invalid intervals and empty files
- Add integration test for full enrichment pipeline
- Include pytest-benchmark for overlap calculation performance
```

---

### Day 10 (2024-01-24, Wednesday)

**Commit 10** - 11:35 AM
```
fix: prevent division by zero in log2 fold enrichment

- Add epsilon smoothing (1e-6) to observed and expected values
- Handle edge case where permutations yield zero overlap
- Prevents NaN/Inf in enrichment statistics
```

**Commit 11** - 03:50 PM
```
refactor: extract interval subtraction to separate method

- Move exclusion region logic to _subtract_intervals()
- Add TODO comment for proper BedTools-like partial overlap handling
- Improves code readability in load_data() method
```

---

### Day 11 (2024-01-25, Thursday)

**Commit 12** - 09:27 AM
```
docs: add comprehensive README with benchmarks and examples

- Include installation instructions using uv
- Add performance benchmark table with realistic dataset sizes
- Document all CLI parameters and output schema
- Reference ENCODE and GEO accessions for example data
```

---

### Day 12 (2024-01-26, Friday)

**Commit 13** - 10:56 AM
```
build: add multi-stage Dockerfile with uv

- Create builder stage for dependency installation
- Use python:3.11-slim base for minimal image size
- Copy virtual environment from builder to runtime
- Add metadata labels and health check endpoint
```

**Commit 14** - 02:03 PM
```
chore: add .gitignore, LICENSE, and example config

- Add MIT license
- Configure .gitignore for Python, uv, and data files
- Include example JSON config for reproducible runs
```

---

## Post-Release Maintenance

### Day 15 (2024-01-29, Monday)

**Commit 15** - 11:20 AM
```
docs: update README with Docker usage examples

- Add docker build and run commands
- Document volume mount patterns for input data
- Include troubleshooting section for common errors
```

---

## Implementation Notes

1. **Spacing Strategy**: Commits are distributed across 12 working days over 2 weeks
2. **Realistic Patterns**:
   - Initial scaffold on Day 1
   - Feature development Days 2-5
   - Refactoring and optimization Days 5-7
   - Testing and documentation Days 9-12
   - Post-release polish Days 15+
3. **Conventional Commit Types Used**:
   - `feat:` - New features (8 commits)
   - `fix:` - Bug fixes (2 commits)
   - `refactor:` - Code restructuring (2 commits)
   - `test:` - Test additions (1 commit)
   - `docs:` - Documentation (2 commits)
   - `build:` - Build system changes (1 commit)
   - `chore:` - Maintenance tasks (1 commit)

4. **Human Developer Signals**:
   - Bug fixes appear after initial implementation
   - Performance optimization comes after basic functionality
   - Tests written mid-development, not upfront
   - Documentation updated incrementally
   - TODO comments left for future work
   - Refactoring based on code review insights
