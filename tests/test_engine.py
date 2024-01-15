"""
Unit tests and benchmarks for rloop-eclip-engine.

Tests cover:
- Chromosome normalization edge cases
- Empty input file handling
- Interval overlap calculation
- Permutation test correctness
- Configuration validation
"""

import tempfile
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from rloop_eclip.core import ChromosomeNormalizer, calculate_total_overlap, load_bed_to_polars
from rloop_eclip.engine import EnrichmentEngine
from rloop_eclip.schemas import BEDRecord, EnrichmentConfig


class TestChromosomeNormalizer:
    """Test suite for chromosome name normalization."""

    def test_ucsc_format_preserved(self):
        """UCSC format names should pass through unchanged."""
        assert ChromosomeNormalizer.normalize_to_ucsc("chr1") == "chr1"
        assert ChromosomeNormalizer.normalize_to_ucsc("chrX") == "chrX"
        assert ChromosomeNormalizer.normalize_to_ucsc("chrM") == "chrM"

    def test_ensembl_to_ucsc_conversion(self):
        """Ensembl names should convert to UCSC format."""
        assert ChromosomeNormalizer.normalize_to_ucsc("1") == "chr1"
        assert ChromosomeNormalizer.normalize_to_ucsc("X") == "chrX"
        assert ChromosomeNormalizer.normalize_to_ucsc("MT") == "chrM"

    def test_ucsc_to_ensembl_conversion(self):
        """UCSC names should convert to Ensembl format."""
        assert ChromosomeNormalizer.normalize_to_ensembl("chr1") == "1"
        assert ChromosomeNormalizer.normalize_to_ensembl("chrX") == "X"
        assert ChromosomeNormalizer.normalize_to_ensembl("chrM") == "MT"

    def test_exclusion_patterns(self):
        """Unlocalized and random contigs should be excluded."""
        assert ChromosomeNormalizer.should_exclude("chr1_random")
        assert ChromosomeNormalizer.should_exclude("chrUn_GL000220v1")
        assert ChromosomeNormalizer.should_exclude("chr1_gl000191_alt")
        assert ChromosomeNormalizer.should_exclude("chrM")
        assert not ChromosomeNormalizer.should_exclude("chr1")


class TestBEDParsing:
    """Test suite for BED file loading and validation."""

    def test_load_valid_bed3(self, tmp_path):
        """Load valid BED3 format file."""
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t1000\t2000\nchr2\t5000\t6000\n")

        lf = load_bed_to_polars(bed_file, expected_format="bed3")
        df = lf.collect()

        assert df.height == 2
        assert df["chrom"].to_list() == ["chr1", "chr2"]
        assert df["start"].to_list() == [1000, 5000]

    def test_invalid_intervals_filtered(self, tmp_path):
        """Intervals with start >= end should be filtered out."""
        bed_file = tmp_path / "invalid.bed"
        bed_file.write_text(
            "chr1\t1000\t2000\n"
            "chr1\t3000\t3000\n"  # start == end (invalid)
            "chr1\t5000\t4000\n"  # start > end (invalid)
        )

        lf = load_bed_to_polars(bed_file, expected_format="bed3")
        df = lf.collect()

        # Only first interval should remain
        assert df.height == 1
        assert df["start"][0] == 1000

    def test_empty_bed_file_handling(self, tmp_path):
        """Empty BED files should result in empty LazyFrame without crash."""
        bed_file = tmp_path / "empty.bed"
        bed_file.write_text("# Comment line only\n")

        lf = load_bed_to_polars(bed_file, expected_format="bed3")
        df = lf.collect()

        assert df.height == 0

    def test_chromosome_normalization_applied(self, tmp_path):
        """Ensembl chromosome names should be normalized to UCSC."""
        bed_file = tmp_path / "ensembl.bed"
        bed_file.write_text("1\t1000\t2000\nX\t5000\t6000\n")

        lf = load_bed_to_polars(bed_file, expected_format="bed3", normalize_chroms=True)
        df = lf.collect()

        assert df["chrom"].to_list() == ["chr1", "chrX"]


class TestIntervalOverlap:
    """Test suite for interval overlap calculations."""

    def test_simple_overlap(self):
        """Test basic overlapping intervals."""
        intervals_a = pl.DataFrame({
            "chrom": ["chr1"],
            "start": [1000],
            "end": [2000],
        }).lazy()

        intervals_b = pl.DataFrame({
            "chrom": ["chr1"],
            "start": [1500],
            "end": [2500],
        }).lazy()

        overlap = calculate_total_overlap(intervals_a, intervals_b)
        assert overlap == 500  # 1500-2000 overlap

    def test_no_overlap_different_chromosomes(self):
        """Intervals on different chromosomes should have zero overlap."""
        intervals_a = pl.DataFrame({
            "chrom": ["chr1"],
            "start": [1000],
            "end": [2000],
        }).lazy()

        intervals_b = pl.DataFrame({
            "chrom": ["chr2"],
            "start": [1000],
            "end": [2000],
        }).lazy()

        overlap = calculate_total_overlap(intervals_a, intervals_b)
        assert overlap == 0

    def test_no_overlap_same_chromosome(self):
        """Non-overlapping intervals on same chromosome should have zero overlap."""
        intervals_a = pl.DataFrame({
            "chrom": ["chr1"],
            "start": [1000],
            "end": [2000],
        }).lazy()

        intervals_b = pl.DataFrame({
            "chrom": ["chr1"],
            "start": [3000],
            "end": [4000],
        }).lazy()

        overlap = calculate_total_overlap(intervals_a, intervals_b)
        assert overlap == 0


class TestEnrichmentEngine:
    """Integration tests for complete enrichment analysis."""

    @pytest.fixture
    def minimal_test_data(self, tmp_path):
        """Create minimal test dataset for integration testing."""
        # Create eCLIP BED file
        eclip_bed = tmp_path / "eclip.bed"
        eclip_bed.write_text(
            "chr1\t10000\t11000\tpeak1\t100\t+\n"
            "chr1\t20000\t21000\tpeak2\t100\t+\n"
            "chr2\t30000\t31000\tpeak3\t100\t+\n"
        )

        # Create DRIPc BED file (overlaps with first eCLIP peak)
        drip_bed = tmp_path / "drip.bed"
        drip_bed.write_text(
            "chr1\t10500\t12000\trloop1\t50\t+\n"
            "chr2\t35000\t36000\trloop2\t50\t+\n"
        )

        # Create genome sizes
        genome_sizes = tmp_path / "genome.sizes"
        genome_sizes.write_text("chr1\t100000\nchr2\t100000\n")

        return {
            "eclip_bed": eclip_bed,
            "drip_bed": drip_bed,
            "genome_sizes": genome_sizes,
        }

    def test_enrichment_engine_execution(self, minimal_test_data, tmp_path):
        """Test complete enrichment analysis with minimal dataset."""
        output = tmp_path / "results.parquet"

        config = EnrichmentConfig(
            eclip_bed=minimal_test_data["eclip_bed"],
            drip_bed=minimal_test_data["drip_bed"],
            genome_sizes=minimal_test_data["genome_sizes"],
            iterations=100,  # Low iteration count for fast testing
            threads=2,
            output=output,
            seed=42,  # Reproducibility
        )

        engine = EnrichmentEngine(config)
        result = engine.run()

        # Validate result structure
        assert result.n_eclip_peaks == 3
        assert result.n_drip_regions == 2
        assert result.observed_overlap_bp == 500  # chr1:10500-11000
        assert 0.0 <= result.empirical_pvalue <= 1.0
        assert result.n_iterations == 100

    def test_config_validation_invalid_file(self, tmp_path):
        """Configuration should fail for non-existent input files."""
        with pytest.raises(ValueError, match="does not exist"):
            EnrichmentConfig(
                eclip_bed=Path("/nonexistent/file.bed"),
                drip_bed=tmp_path / "drip.bed",
                genome_sizes=tmp_path / "genome.sizes",
                output=tmp_path / "out.parquet",
            )

    def test_config_validation_invalid_output_format(self, tmp_path):
        """Output must be .parquet or .json format."""
        eclip_bed = tmp_path / "eclip.bed"
        eclip_bed.touch()

        with pytest.raises(ValueError, match="must be .parquet or .json"):
            EnrichmentConfig(
                eclip_bed=eclip_bed,
                drip_bed=eclip_bed,
                genome_sizes=tmp_path / "genome.sizes",
                output=tmp_path / "output.txt",  # Invalid format
            )


class TestBEDRecordSchema:
    """Test Pydantic BED record validation."""

    def test_valid_bed_record(self):
        """Valid BED record should pass validation."""
        record = BEDRecord(chrom="chr1", start=1000, end=2000, strand="+")
        assert record.chrom == "chr1"
        assert record.start == 1000
        assert record.end == 2000

    def test_invalid_interval_rejected(self):
        """Intervals with start >= end should be rejected."""
        with pytest.raises(ValueError, match="Invalid interval"):
            BEDRecord(chrom="chr1", start=2000, end=1000)

        with pytest.raises(ValueError, match="Invalid interval"):
            BEDRecord(chrom="chr1", start=1000, end=1000)


# Benchmark tests (run with: pytest --benchmark-only)
@pytest.mark.benchmark
def test_overlap_calculation_benchmark(benchmark):
    """Benchmark interval overlap calculation performance."""
    # Generate realistic dataset size
    n_intervals = 10000

    intervals_a = pl.DataFrame({
        "chrom": [f"chr{i % 22 + 1}" for i in range(n_intervals)],
        "start": np.random.randint(0, 100000000, n_intervals),
        "end": np.random.randint(100000000, 200000000, n_intervals),
    }).lazy()

    intervals_b = pl.DataFrame({
        "chrom": [f"chr{i % 22 + 1}" for i in range(n_intervals)],
        "start": np.random.randint(0, 100000000, n_intervals),
        "end": np.random.randint(100000000, 200000000, n_intervals),
    }).lazy()

    result = benchmark(calculate_total_overlap, intervals_a, intervals_b)
    assert result >= 0  # Basic sanity check
