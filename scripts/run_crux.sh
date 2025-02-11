#!/bin/bash
set -e

echo "[INFO]START: $(date)"
WORK_DIR=$(realpath "$PWD/..")
FASTA="${WORK_DIR}/data/target.fasta"
MZML="${WORK_DIR}/data/PDL-2-5.mzML"
OUTPUT_DIR="${WORK_DIR}/output/crux"

if [ ! -f "$FASTA" ]; then
  echo "Error: FASTA file not found at $FASTA"
  exit 1
fi

if [ ! -f "$MZML" ]; then
  echo "Error: mzML file not found at $MZML"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Run Crux Tide Index
INDEX_DIR="$OUTPUT_DIR/tide-index"
echo "[INFO]Running crux tide-index..."
mkdir -p "$INDEX_DIR"
crux tide-index \
  --enzyme trypsin \
  --missed-cleavages 2 \
  --decoy-format peptide-reverse \
  --peptide-list T \
  --output-dir "$INDEX_DIR" \
  --overwrite T \
  --mods-spec 4M+15.9949 \
  "$FASTA" "$OUTPUT_DIR/idx"

# Run Crux Tide Search
SEARCH_DIR="$OUTPUT_DIR/tide-search"
echo "[INFO]Running crux tide-search..."
mkdir -p "$SEARCH_DIR"
crux tide-search \
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
PERC_DIR="$OUTPUT_DIR/crux_perc"
echo "[INFO]Running crux percolator..."
mkdir -p "$PERC_DIR"
crux percolator \
  --top-match 1 \
  --output-dir "$PERC_DIR" \
  --overwrite T \
  "$SEARCH_DIR/tide-search.txt"

# Copy the pin tsv file to the data directory
echo "[INFO]Copying pin tsv file to data directory..."
cp "$PERC_DIR/make-pin.pin" "$WORK_DIR/data/"

echo "[INFO]END: $(date)"