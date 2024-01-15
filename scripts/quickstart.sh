#!/bin/bash
# Quick start script for rloop-eclip-engine demonstration
# This script assumes you have sample data in the data/ directory

set -e

echo "================================"
echo "rloop-eclip-engine Quick Start"
echo "================================"
echo ""

# Check if virtual environment is activated
if [ -z "$VIRTUAL_ENV" ]; then
    echo "Error: Virtual environment not activated"
    echo "Run: source .venv/bin/activate"
    exit 1
fi

# Check for required input files
ECLIP_BED="data/sample_eclip.bed"
DRIP_BED="data/sample_drip.bed"
GENOME_SIZES="data/hg38.chrom.sizes"

if [ ! -f "$ECLIP_BED" ]; then
    echo "Error: eCLIP BED file not found: $ECLIP_BED"
    echo "Please download sample data first"
    exit 1
fi

if [ ! -f "$DRIP_BED" ]; then
    echo "Error: DRIPc BED file not found: $DRIP_BED"
    echo "Please download sample data first"
    exit 1
fi

if [ ! -f "$GENOME_SIZES" ]; then
    echo "Error: Genome sizes file not found: $GENOME_SIZES"
    echo "Download from: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
    exit 1
fi

# Create output directory
mkdir -p results

# Run analysis
echo "Running enrichment analysis..."
echo "  eCLIP peaks: $ECLIP_BED"
echo "  DRIPc regions: $DRIP_BED"
echo "  Iterations: 1000 (quick test)"
echo ""

rloop-eclip \
    --eclip-bed "$ECLIP_BED" \
    --drip-bed "$DRIP_BED" \
    --genome-sizes "$GENOME_SIZES" \
    --iterations 1000 \
    --threads 4 \
    --output results/quickstart_results.parquet \
    --seed 42

echo ""
echo "================================"
echo "Analysis complete!"
echo "Results saved to: results/quickstart_results.parquet"
echo "================================"
