#!/bin/bash
set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 OUT_DIR SEED"
    exit 1
fi

OUT_DIR=$1
SEED=$2

echo "[INFO]START: $(date)"
DATA_DIR="/proj/proteoforma_nsc/smooth_q_to_pep/datasets/PXD004948"
FASTA="${DATA_DIR}/Mus-musculus.fasta"
MZML="${DATA_DIR}/20160323_CoAN_CTRL1_3372.mzML"
WORK_DIR="/proj/proteoforma_nsc/smooth_q_to_pep/${OUT_DIR}/PXD004948/run${SEED}"
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
  --min-length 7 \
  --decoy-format shuffle \
  --peptide-list T \
  --output-dir "$INDEX_DIR" \
  --overwrite T \
  --max-mods 3 \
  --mods-spec 1M+15.994915,1N+0.984016,1Q+0.984016 \
  --nterm-peptide-mods-spec 1X+42.010565,1X+43.005814,1X-17.026549,1X+25.980265 \
  --seed $SEED \
  "$FASTA" "$OUTPUT_DIR/idx"

# Run Crux Tide Search
SEARCH_DIR="$OUTPUT_DIR/tide-search"
echo "[INFO]Running crux tide-search..."
mkdir -p "$SEARCH_DIR"
/proj/proteoforma_nsc/crux tide-search \
  --top-match 1 \
  --precursor-window 4.5 \
  --precursor-window-type ppm \
  --fragment-tolerance 0.02 \
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
