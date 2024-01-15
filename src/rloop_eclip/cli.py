"""
Command-line interface for rloop-eclip-engine.

Uses Click for argument parsing and Rich for terminal formatting.
"""

import sys
from pathlib import Path

import click
from loguru import logger
from pydantic import ValidationError
from rich.console import Console
from rich.table import Table

from .engine import EnrichmentEngine
from .schemas import EnrichmentConfig, EnrichmentResult


# Configure loguru to output structured logs
logger.remove()  # Remove default handler
logger.add(
    sys.stderr,
    format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
    level="INFO",
)

console = Console()


@click.command()
@click.option(
    "--eclip-bed",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to eCLIP peak BED file (BED6 or narrowPeak format)",
)
@click.option(
    "--drip-bed",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to DRIPc-seq R-loop BED file",
)
@click.option(
    "--genome-sizes",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Chromosome sizes file (e.g., hg38.chrom.sizes)",
)
@click.option(
    "--exclude-regions",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Optional ENCODE blacklist BED file",
)
@click.option(
    "--iterations",
    type=int,
    default=10000,
    show_default=True,
    help="Number of Monte Carlo permutations",
)
@click.option(
    "--threads",
    type=int,
    default=4,
    show_default=True,
    help="Number of CPU threads for Polars parallelism",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output file path (.parquet or .json)",
)
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility",
)
@click.option(
    "--verbose",
    is_flag=True,
    help="Enable verbose logging (DEBUG level)",
)
def main(
    eclip_bed: Path,
    drip_bed: Path,
    genome_sizes: Path,
    exclude_regions: Path,
    iterations: int,
    threads: int,
    output: Path,
    seed: int,
    verbose: bool,
):
    """
    rloop-eclip-engine: High-performance genomic interval enrichment analysis.

    Quantifies spatial enrichment between eCLIP RBP peaks and DRIPc-seq R-loop regions
    using vectorized Monte Carlo permutation testing.

    Example:

        rloop-eclip --eclip-bed DHX9_peaks.bed --drip-bed rloops.bed \\
                    --genome-sizes hg38.chrom.sizes --iterations 10000 \\
                    --output results.parquet
    """
    # Enable verbose logging if requested
    if verbose:
        logger.remove()
        logger.add(sys.stderr, level="DEBUG")
        logger.debug("Verbose logging enabled")

    console.print("[bold cyan]rloop-eclip-engine v0.2.1[/bold cyan]")
    console.print("High-performance genomic interval enrichment analysis\n")

    try:
        # Validate configuration using Pydantic
        config = EnrichmentConfig(
            eclip_bed=eclip_bed,
            drip_bed=drip_bed,
            genome_sizes=genome_sizes,
            exclude_regions=exclude_regions,
            iterations=iterations,
            threads=threads,
            output=output,
            seed=seed,
        )

        # Initialize engine and run analysis
        engine = EnrichmentEngine(config)
        result = engine.run()

        # Save results
        save_results(result, output)

        # Display formatted results table
        display_results(result)

        console.print(f"\n[bold green]Results saved to:[/bold green] {output}")

    except ValidationError as e:
        console.print(f"[bold red]Configuration validation failed:[/bold red]")
        for error in e.errors():
            field = " -> ".join(str(loc) for loc in error["loc"])
            console.print(f"  • {field}: {error['msg']}")
        sys.exit(1)

    except Exception as e:
        logger.exception("Fatal error during execution")
        console.print(f"[bold red]Error:[/bold red] {e}")
        sys.exit(1)


def save_results(result: EnrichmentResult, output_path: Path) -> None:
    """Save enrichment results to Parquet or JSON format."""
    import json

    result_dict = result.model_dump()

    if output_path.suffix == ".parquet":
        # Save as Parquet using Polars for consistency
        import polars as pl

        df = pl.DataFrame([result_dict])
        df.write_parquet(output_path)
        logger.info(f"Results written to Parquet: {output_path}")

    elif output_path.suffix == ".json":
        with open(output_path, "w") as f:
            json.dump(result_dict, f, indent=2)
        logger.info(f"Results written to JSON: {output_path}")

    else:
        raise ValueError(f"Unsupported output format: {output_path.suffix}")


def display_results(result: EnrichmentResult) -> None:
    """Display formatted results table in terminal using Rich."""
    table = Table(title="Enrichment Analysis Results", show_header=True, header_style="bold magenta")
    table.add_column("Metric", style="cyan", width=30)
    table.add_column("Value", style="green", justify="right")

    # Format results with appropriate precision
    table.add_row("Observed Overlap (bp)", f"{result.observed_overlap_bp:,}")
    table.add_row("Expected Overlap (bp)", f"{result.expected_overlap_bp:,.1f}")
    table.add_row("Log2 Fold Enrichment", f"{result.log2_fold_enrichment:.3f}")
    table.add_row("Empirical p-value", f"{result.empirical_pvalue:.4e}")
    table.add_row("95% CI Lower", f"{result.ci_95_lower:,.1f}")
    table.add_row("95% CI Upper", f"{result.ci_95_upper:,.1f}")
    table.add_row("eCLIP Peaks", f"{result.n_eclip_peaks:,}")
    table.add_row("R-loop Regions", f"{result.n_drip_regions:,}")
    table.add_row("Permutations", f"{result.n_iterations:,}")
    table.add_row("Runtime (seconds)", f"{result.runtime_seconds:.2f}")

    console.print(table)


if __name__ == "__main__":
    main()
