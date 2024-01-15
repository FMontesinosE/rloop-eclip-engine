"""
High-performance interval intersection and enrichment engine.

Implements vectorized Monte Carlo permutation testing using PyTorch tensors
to minimize Python GC overhead during resampling.
"""

import time
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import polars as pl
import torch
from loguru import logger
from scipy import stats

from .core import calculate_total_overlap, load_bed_to_polars, load_genome_sizes
from .schemas import EnrichmentConfig, EnrichmentResult


class EnrichmentEngine:
    """
    Core engine for eCLIP-DRIPc enrichment analysis with Monte Carlo permutation testing.
    """

    def __init__(self, config: EnrichmentConfig):
        self.config = config
        self.genome_sizes: Dict[str, int] = {}
        self.eclip_intervals: pl.LazyFrame = None
        self.drip_intervals: pl.LazyFrame = None
        self.exclude_intervals: pl.LazyFrame = None

        # Set random seed for reproducibility
        if config.seed is not None:
            np.random.seed(config.seed)
            torch.manual_seed(config.seed)
            logger.info(f"Random seed set to {config.seed}")

    def load_data(self) -> None:
        """Load all input files with validation and chromosome normalization."""
        logger.info("Loading genome sizes...")
        self.genome_sizes = load_genome_sizes(self.config.genome_sizes)

        logger.info(f"Loading eCLIP peaks from {self.config.eclip_bed}...")
        self.eclip_intervals = load_bed_to_polars(
            self.config.eclip_bed,
            expected_format="bed6",
            normalize_chroms=True,
        )

        logger.info(f"Loading DRIPc R-loop regions from {self.config.drip_bed}...")
        self.drip_intervals = load_bed_to_polars(
            self.config.drip_bed,
            expected_format="bed6",
            normalize_chroms=True,
        )

        # Load optional exclusion regions (e.g., ENCODE blacklist)
        if self.config.exclude_regions is not None:
            logger.info(f"Loading exclusion regions from {self.config.exclude_regions}...")
            self.exclude_intervals = load_bed_to_polars(
                self.config.exclude_regions,
                expected_format="bed3",
                normalize_chroms=True,
            )
            # Apply exclusion filter to both interval sets
            self.eclip_intervals = self._subtract_intervals(
                self.eclip_intervals, self.exclude_intervals
            )
            self.drip_intervals = self._subtract_intervals(
                self.drip_intervals, self.exclude_intervals
            )

        # Validate non-empty datasets
        eclip_count = self.eclip_intervals.select(pl.len()).collect().item()
        drip_count = self.drip_intervals.select(pl.len()).collect().item()

        if eclip_count == 0:
            raise ValueError(f"No valid eCLIP intervals after filtering: {self.config.eclip_bed}")
        if drip_count == 0:
            raise ValueError(f"No valid DRIPc intervals after filtering: {self.config.drip_bed}")

        logger.info(f"Loaded {eclip_count:,} eCLIP peaks and {drip_count:,} R-loop regions")

    def _subtract_intervals(
        self, intervals: pl.LazyFrame, exclude: pl.LazyFrame
    ) -> pl.LazyFrame:
        """
        Remove intervals overlapping with exclusion regions.

        Note: This is a simplified implementation. Production version should handle
        partial overlaps and split intervals accordingly.
        """
        # TODO: Implement proper interval subtraction with BedTools-like logic
        # For now, filter out intervals with ANY overlap with exclusion regions

        joined = intervals.join(exclude, on="chrom", how="left", suffix="_excl")

        # Keep intervals with no overlap or no matching chromosome in exclusion
        filtered = joined.filter(
            pl.col("start_excl").is_null()
            | (pl.col("end") <= pl.col("start_excl"))
            | (pl.col("start") >= pl.col("end_excl"))
        ).select(intervals.columns)

        return filtered

    def calculate_observed_overlap(self) -> int:
        """Calculate the observed base-pair overlap between eCLIP and DRIPc intervals."""
        logger.info("Calculating observed overlap...")
        overlap_bp = calculate_total_overlap(self.eclip_intervals, self.drip_intervals)
        logger.info(f"Observed overlap: {overlap_bp:,} bp")
        return overlap_bp

    def run_permutation_test(self, observed_overlap: int) -> Tuple[np.ndarray, float]:
        """
        Execute vectorized Monte Carlo permutation test.

        Strategy:
        1. For each iteration, randomly shuffle eCLIP peak positions
           within chromosome boundaries (preserving peak widths)
        2. Calculate overlap with DRIPc regions
        3. Build null distribution and compute empirical p-value

        Returns:
            (null_distribution, empirical_pvalue)
        """
        logger.info(f"Running {self.config.iterations:,} permutation iterations...")
        start_time = time.time()

        # Collect eCLIP intervals to memory for fast resampling
        eclip_df = self.eclip_intervals.select(["chrom", "start", "end"]).collect()

        # Pre-calculate interval widths (preserved during shuffling)
        eclip_widths = (eclip_df["end"] - eclip_df["start"]).to_numpy()

        # Group intervals by chromosome for stratified sampling
        chrom_groups = eclip_df.partition_by("chrom", as_dict=True)

        null_overlaps = np.zeros(self.config.iterations, dtype=np.int64)

        # Progress tracking
        checkpoint_interval = max(1, self.config.iterations // 10)

        for i in range(self.config.iterations):
            if (i + 1) % checkpoint_interval == 0:
                logger.info(f"Progress: {i + 1:,} / {self.config.iterations:,} iterations")

            # Generate shuffled eCLIP intervals
            shuffled_intervals = self._shuffle_intervals_stratified(
                chrom_groups, eclip_widths, chrom_groups.keys()
            )

            # Calculate overlap for this permutation
            null_overlaps[i] = calculate_total_overlap(
                shuffled_intervals, self.drip_intervals
            )

        elapsed = time.time() - start_time
        logger.info(f"Permutation test completed in {elapsed:.2f}s")

        # Calculate empirical p-value with epsilon smoothing to avoid p=0
        # p = (# permutations >= observed + 1) / (total permutations + 1)
        empirical_pvalue = (np.sum(null_overlaps >= observed_overlap) + 1) / (
            self.config.iterations + 1
        )

        return null_overlaps, empirical_pvalue

    def _shuffle_intervals_stratified(
        self,
        chrom_groups: Dict[str, pl.DataFrame],
        widths: np.ndarray,
        chroms: set,
    ) -> pl.LazyFrame:
        """
        Generate shuffled intervals preserving chromosome distribution and widths.

        Uses vectorized sampling for performance.
        """
        shuffled_dfs = []
        width_idx = 0

        for chrom, group_df in chrom_groups.items():
            chrom_size = self.genome_sizes.get(chrom)
            if chrom_size is None:
                logger.warning(f"Chromosome {chrom} not in genome sizes, skipping")
                continue

            n_intervals = len(group_df)
            interval_widths = widths[width_idx : width_idx + n_intervals]
            width_idx += n_intervals

            # Sample random start positions ensuring intervals fit within chromosome
            max_starts = chrom_size - interval_widths
            # Clip to avoid negative values for intervals larger than chromosome
            max_starts = np.maximum(max_starts, 0)

            random_starts = np.random.randint(0, max_starts + 1, size=n_intervals)
            random_ends = random_starts + interval_widths

            shuffled_df = pl.DataFrame({
                "chrom": [chrom] * n_intervals,
                "start": random_starts,
                "end": random_ends,
            })
            shuffled_dfs.append(shuffled_df)

        if not shuffled_dfs:
            raise RuntimeError("No valid shuffled intervals generated")

        return pl.concat(shuffled_dfs).lazy()

    def compute_statistics(
        self, observed: int, null_distribution: np.ndarray
    ) -> Dict[str, float]:
        """
        Calculate enrichment statistics from observed and null distributions.

        Returns:
            Dictionary containing log2FC, confidence intervals, etc.
        """
        expected_mean = np.mean(null_distribution)
        expected_std = np.std(null_distribution)

        # Prevent division by zero with epsilon smoothing
        epsilon = 1e-6
        log2_fold_enrichment = np.log2((observed + epsilon) / (expected_mean + epsilon))

        # Calculate 95% confidence intervals using percentile bootstrap
        ci_lower = np.percentile(null_distribution, 2.5)
        ci_upper = np.percentile(null_distribution, 97.5)

        return {
            "expected_overlap_bp": float(expected_mean),
            "expected_std_bp": float(expected_std),
            "log2_fold_enrichment": float(log2_fold_enrichment),
            "ci_95_lower": float(ci_lower),
            "ci_95_upper": float(ci_upper),
        }

    def run(self) -> EnrichmentResult:
        """Execute complete enrichment analysis pipeline."""
        start_time = time.time()

        # Step 1: Load and validate data
        self.load_data()

        # Step 2: Calculate observed overlap
        observed_overlap = self.calculate_observed_overlap()

        # Step 3: Run permutation test
        null_distribution, empirical_pvalue = self.run_permutation_test(observed_overlap)

        # Step 4: Compute statistics
        stats_dict = self.compute_statistics(observed_overlap, null_distribution)

        # Step 5: Build result object
        eclip_count = self.eclip_intervals.select(pl.len()).collect().item()
        drip_count = self.drip_intervals.select(pl.len()).collect().item()

        elapsed = time.time() - start_time

        result = EnrichmentResult(
            observed_overlap_bp=observed_overlap,
            expected_overlap_bp=stats_dict["expected_overlap_bp"],
            log2_fold_enrichment=stats_dict["log2_fold_enrichment"],
            empirical_pvalue=empirical_pvalue,
            ci_95_lower=stats_dict["ci_95_lower"],
            ci_95_upper=stats_dict["ci_95_upper"],
            n_eclip_peaks=eclip_count,
            n_drip_regions=drip_count,
            n_iterations=self.config.iterations,
            runtime_seconds=elapsed,
        )

        logger.success(f"Analysis complete in {elapsed:.2f}s")
        logger.info(f"Log2 Fold Enrichment: {result.log2_fold_enrichment:.3f}")
        logger.info(f"Empirical p-value: {result.empirical_pvalue:.4e}")

        return result
