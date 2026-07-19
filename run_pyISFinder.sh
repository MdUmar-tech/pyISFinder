#!/bin/bash
# ==============================================================================
# pyISFinder Execution Script for Falsiroseomonas sp. RANM36
# ==============================================================================

set -e

# Work directory
WORKDIR="/Users/mdumar/Desktop/freelancer/Royston/ranm36/draft/ISCSCAN"
cd "$WORKDIR"

echo "======================================================================"
echo "          Running pyISFinder Pipeline"
echo "======================================================================"

# 1. Activate virtual environment
if [ -d "./venv" ]; then
    echo "[1/3] Activating virtual environment..."
    source ./venv/bin/activate
else
    echo "[1/3] Creating virtual environment and installing dependencies..."
    python3 -m venv ./venv
    source ./venv/bin/activate
    pip install -r requirements.txt
fi

# 2. Setup reference database
echo "[2/3] Checking reference database..."
python3 download_db.py

# 3. Execute pyISFinder pipeline on genome
echo "[3/3] Running pyISFinder analysis on Falsiroseomonas_sp._RANM36.fasta..."
python3 pyISFinder.py -i Falsiroseomonas_sp._RANM36.fasta -o results

echo ""
echo "======================================================================"
echo "SUCCESS! Results generated in: $WORKDIR/results"
echo "  - Predictions TSV: results/Falsiroseomonas_sp._RANM36_IS_predictions.tsv"
echo "  - GFF3 Track:      results/Falsiroseomonas_sp._RANM36_IS.gff3"
echo "  - Summary Table:   results/Falsiroseomonas_sp._RANM36_IS_summary.txt"
echo "  - Genomic Map:     results/Falsiroseomonas_sp._RANM36_IS_genome_map.png"
echo "======================================================================"
