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
ISOTONIC_DIR="$WORK_DIR/$OUTPUT_DIR/IsoLogReg"
ISOTONIC_NORMAL_DIR="$ISOTONIC_DIR/0"
IP_ISOTONIC_DIR="$WORK_DIR/$OUTPUT_DIR/ipPEP"
IP_ISOTONIC_NORMAL_DIR="$IP_ISOTONIC_DIR/0"
DATA_OUTPUT_DIR="/data/$OUTPUT_DIR"

cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }
# Create necessary directories
mkdir -p "$OUTPUT_DIR" "$ISOTONIC_DIR" "$ISOTONIC_NORMAL_DIR" "$IP_ISOTONIC_DIR" "$IP_ISOTONIC_NORMAL_DIR"

# Run isotonic Percolator with isotonic logistic regression
echo "[INFO]Running Percolator with isotonic logistic regression..."
podman run --rm -v "$WORK_DIR:/data" percolator:isotonic \
    /usr/bin/percolator -S $SEED -Y \
    -w "$DATA_OUTPUT_DIR/IsoLogReg/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/IsoLogReg/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/IsoLogReg/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/IsoLogReg/0/peptide.decoy.txt" \
    /data/data/make-pin.pin 2>&1 | tee "$ISOTONIC_NORMAL_DIR/percolator.log"

# Run isotonic Percolator with isotonic interpolation
echo "[INFO]Running Percolator with isotonic interpolation..."
podman run --rm -v "$WORK_DIR:/data" percolator:isotonic \
    /usr/bin/percolator -S $SEED --ip-pep -Y \
    -w "$DATA_OUTPUT_DIR/ipPEP/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/ipPEP/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/ipPEP/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/ipPEP/0/peptide.decoy.txt" \
    /data/data/make-pin.pin 2>&1 | tee "$IP_ISOTONIC_NORMAL_DIR/percolator.log"

echo "[INFO]END: $(date)"