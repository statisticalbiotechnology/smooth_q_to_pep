#!/bin/bash
set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 OUT_DIR SEED"
    exit 1
fi

OUT_DIR=$1
SEED=$2

echo "[INFO]START: $(date)"
DATA_DIR="/proj/proteoforma_nsc/smooth_q_to_pep/datasets/PXD013274"
FASTA="${DATA_DIR}/target.fasta"
MZML="${DATA_DIR}/PDL-2-5.mzML"
WORK_DIR="/proj/proteoforma_nsc/smooth_q_to_pep/${OUT_DIR}/PXD013274/run${SEED}"
OUTPUT_DIR="${WORK_DIR}/crux"
INTERM_DIR="${WORK_DIR}/interm"

if [ ! -f "$FASTA" ]; then
  echo "Error: FASTA file not found at $FASTA"
  exit 1
fi

if [ ! -f "$MZML" ]; then
  echo "Error: mzML file not found at $MZML"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$INTERM_DIR"

# Run Crux Tide Index
INDEX_DIR="$OUTPUT_DIR/tide-index"
echo "[INFO]Running crux tide-index..."
mkdir -p "$INDEX_DIR"
/proj/proteoforma_nsc/crux tide-index \
  --enzyme trypsin \
  --missed-cleavages 2 \
  --decoy-format shuffle \
  --peptide-list T \
  --output-dir "$INDEX_DIR" \
  --overwrite T \
  --mods-spec 4M+15.9949 \
  --seed $SEED \
  "$FASTA" "$OUTPUT_DIR/idx"

# Run Crux Tide Search
SEARCH_DIR="$OUTPUT_DIR/tide-search"
echo "[INFO]Running crux tide-search..."
mkdir -p "$SEARCH_DIR"
/proj/proteoforma_nsc/crux tide-search \
  --top-match 1 \
  --precursor-window 10.0 \
  --precursor-window-type ppm \
  --fragment-tolerance 0.02 \
  --min-precursor-charge 2 \
  --use-flanking-peaks true \
  --output-dir "$SEARCH_DIR" \
  --overwrite T \
  --concat T \
  "$MZML" "$OUTPUT_DIR/idx"

# Run Crux Percolator
PERC_DIR="$OUTPUT_DIR/percolator"
echo "[INFO]Running crux percolator..."
mkdir -p "$PERC_DIR"
/proj/proteoforma_nsc/crux percolator \
  --top-match 1 \
  --output-dir "$PERC_DIR" \
  --overwrite T \
  "$SEARCH_DIR/tide-search.txt"

# Copy the pin tsv file to the interm directory
echo "[INFO]Copying pin tsv file to interm data directory..."
cp "$PERC_DIR/make-pin.pin" "$INTERM_DIR/"

echo "[INFO]END: $(date)"