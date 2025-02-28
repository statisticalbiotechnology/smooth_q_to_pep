#!/bin/bash
set -e

echo "[INFO]START: $(date)"
SEED=1
OUTPUT_DIR="output"
if [ ! -z "$1" ]; then
    SEED=$1
fi
if [ ! -z "$2" ]; then
    OUTPUT_DIR=$2
fi

WORK_DIR=$(realpath "$PWD/..")

# Define base directories
SCRIPTS_DIR="$WORK_DIR/scripts"
ISOTONIC_DIR="$WORK_DIR/$OUTPUT_DIR/isoPEP"
ISOTONIC_NORMAL_DIR="$ISOTONIC_DIR/0"
IP_ISOTONIC_DIR="$WORK_DIR/$OUTPUT_DIR/ipPEP"
IP_ISOTONIC_NORMAL_DIR="$IP_ISOTONIC_DIR/0"
TDC_DIR="$WORK_DIR/$OUTPUT_DIR/tdcPEP"
TDC_NORMAL_DIR="$TDC_DIR/0"
DATA_OUTPUT_DIR="/data/$OUTPUT_DIR"

cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }
# Create necessary directories
mkdir -p "$OUTPUT_DIR" "$ISOTONIC_DIR" "$ISOTONIC_NORMAL_DIR" "$IP_ISOTONIC_DIR" "$IP_ISOTONIC_NORMAL_DIR" "$TDC_DIR" "$TDC_NORMAL_DIR"

echo "[INFO]Running isotonic Percolator using TDC..."
podman run --rm -v "$WORK_DIR:/data" percolator:master \
    /usr/bin/percolator -S $SEED -Y \
    -w "$DATA_OUTPUT_DIR/tdcPEP/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/tdcPEP/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/tdcPEP/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/tdcPEP/0/peptide.decoy.txt" \
    /data/data/make-pin.pin 2>&1 | tee "$TDC_NORMAL_DIR/percolator.log"

echo "[INFO]Running isotonic Percolator using isotonic regression..."
podman run --rm -v "$WORK_DIR:/data" percolator:master \
    /usr/bin/percolator -S $SEED --pep-from-q -Y \
    -w "$DATA_OUTPUT_DIR/isoPEP/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/isoPEP/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/isoPEP/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/isoPEP/0/peptide.decoy.txt" \
    /data/data/make-pin.pin 2>&1 | tee "$ISOTONIC_NORMAL_DIR/percolator.log"

# Run isotonic Percolator with isotonic interpolation
echo "[INFO]Running isotonic Percolator using isotonic regression with interpolation...."
podman run --rm -v "$WORK_DIR:/data" percolator:master \
    /usr/bin/percolator -S $SEED --ip-pep -Y \
    -w "$DATA_OUTPUT_DIR/ipPEP/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/ipPEP/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/ipPEP/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/ipPEP/0/peptide.decoy.txt" \
    /data/data/make-pin.pin 2>&1 | tee "$IP_ISOTONIC_NORMAL_DIR/percolator.log"

echo "[INFO]END: $(date)"