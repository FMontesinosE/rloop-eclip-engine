"""
Pydantic v2 schemas for strict input validation and configuration management.
"""

from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, Field, field_validator, model_validator


class EnrichmentConfig(BaseModel):
    """Configuration schema for enrichment analysis execution."""

    eclip_bed: Path = Field(..., description="Path to eCLIP peak BED file")
    drip_bed: Path = Field(..., description="Path to DRIPc-seq R-loop BED file")
    genome_sizes: Path = Field(..., description="Chromosome sizes file")
    exclude_regions: Optional[Path] = Field(None, description="Optional blacklist BED file")
    iterations: int = Field(10000, ge=100, le=1000000, description="Monte Carlo iterations")
    threads: int = Field(4, ge=1, le=128, description="Number of CPU threads")
    output: Path = Field(..., description="Output file path (.parquet or .json)")
    seed: Optional[int] = Field(None, description="Random seed for reproducibility")

    @field_validator("eclip_bed", "drip_bed", "genome_sizes")
    @classmethod
    def validate_input_exists(cls, v: Path) -> Path:
        """Ensure input files exist before execution."""
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        if not v.is_file():
            raise ValueError(f"Path is not a file: {v}")
        return v

    @field_validator("exclude_regions")
    @classmethod
    def validate_optional_file(cls, v: Optional[Path]) -> Optional[Path]:
        """Validate optional exclusion BED file if provided."""
        if v is not None and not v.exists():
            raise ValueError(f"Exclusion file does not exist: {v}")
        return v

    @field_validator("output")
    @classmethod
    def validate_output_format(cls, v: Path) -> Path:
        """Ensure output has valid extension."""
        valid_extensions = {".parquet", ".json"}
        if v.suffix not in valid_extensions:
            raise ValueError(f"Output must be .parquet or .json, got: {v.suffix}")
        return v

    @model_validator(mode="after")
    def validate_thread_count(self) -> "EnrichmentConfig":
        """Warn if thread count exceeds iterations (inefficient)."""
        if self.threads > self.iterations:
            # Log warning but don't fail - user may have specific reasons
            pass
        return self


class EnrichmentResult(BaseModel):
    """Schema for structured enrichment analysis results."""

    observed_overlap_bp: int = Field(..., description="Observed overlap in base pairs")
    expected_overlap_bp: float = Field(..., description="Mean expected overlap from permutations")
    log2_fold_enrichment: float = Field(..., description="Log2(observed / expected)")
    empirical_pvalue: float = Field(..., ge=0.0, le=1.0, description="Empirical p-value")
    ci_95_lower: float = Field(..., description="95% CI lower bound")
    ci_95_upper: float = Field(..., description="95% CI upper bound")
    n_eclip_peaks: int = Field(..., description="Number of eCLIP intervals")
    n_drip_regions: int = Field(..., description="Number of R-loop regions")
    n_iterations: int = Field(..., description="Permutations executed")
    runtime_seconds: float = Field(..., description="Total execution time")

    model_config = {"frozen": True}  # Immutable results


class BEDRecord(BaseModel):
    """Minimal BED format validator for genomic intervals."""

    chrom: str = Field(..., min_length=1, max_length=255)
    start: int = Field(..., ge=0)
    end: int = Field(..., ge=0)
    name: Optional[str] = None
    score: Optional[float] = None
    strand: Optional[Literal["+", "-", "."]] = None

    @model_validator(mode="after")
    def validate_interval(self) -> "BEDRecord":
        """Ensure start < end for valid genomic interval."""
        if self.start >= self.end:
            raise ValueError(f"Invalid interval: start={self.start} >= end={self.end}")
        return self
