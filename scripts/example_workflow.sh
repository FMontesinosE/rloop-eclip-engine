#!/bin/bash
# Example workflow for eCLIP vs DRIPc-seq enrichment analysis
# Demonstrates complete analysis pipeline from data QC to interpretation

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "================================================================"
echo " rloop-eclip-engine: Complete Analysis Workflow"
echo "================================================================"
echo ""

# Check if virtual environment is activated
if [ -z "$VIRTUAL_ENV" ]; then
    echo "Activating virtual environment..."
    source .venv/bin/activate
fi

# Step 1: Data Quality Control
echo -e "${BLUE}Step 1: Input Data Summary${NC}"
echo "----------------------------------------------------------------"
python scripts/summarize_bed.py \
    data/DICER_eCLIP_peaks.bed \
    data/DMSO_DRIPc_example1.narrowPeak

# Step 2: Run enrichment analysis
echo ""
echo -e "${BLUE}Step 2: Running Enrichment Analysis${NC}"
echo "----------------------------------------------------------------"
echo "Analyzing DICER eCLIP peaks vs DRIPc-seq R-loops..."
echo ""

rloop-eclip \
    --eclip-bed data/DICER_eCLIP_peaks.bed \
    --drip-bed data/DMSO_DRIPc_example1.narrowPeak \
    --genome-sizes data/hg38.chrom.sizes \
    --iterations 1000 \
    --threads 4 \
    --output results/DICER_DRIPc_example1_1k.parquet \
    --seed 42

# Step 3: Interpret results
echo ""
echo -e "${BLUE}Step 3: Results Interpretation${NC}"
echo "----------------------------------------------------------------"
python scripts/analyze_results.py results/DICER_DRIPc_example1_1k.parquet

# Step 4: Summary
echo ""
echo -e "${GREEN}✓ Analysis complete!${NC}"
echo ""
echo "Generated files:"
echo "  - results/DICER_DRIPc_example1_1k.parquet   (enrichment statistics)"
echo ""
echo "Next steps:"
echo "  1. Review the empirical p-value and fold enrichment"
echo "  2. Compare with biological controls (if available)"
echo "  3. Visualize overlaps using genome browser"
echo ""
echo "================================================================"
