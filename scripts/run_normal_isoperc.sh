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
ISOTONIC_DIR="$WORK_DIR/$OUTPUT_DIR/isotonic"
ISOTONIC_NORMAL_DIR="$ISOTONIC_DIR/0"
DATA_OUTPUT_DIR="/data/$OUTPUT_DIR"

cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }
# Create necessary directories
mkdir -p "$OUTPUT_DIR" "$ISOTONIC_DIR" "$ISOTONIC_NORMAL_DIR"

# Run isotonic Percolator normally
echo "[INFO]Running isotonic Percolator normally..."
podman run --rm -v "$WORK_DIR:/data" percolator:isotonic \
    /usr/bin/percolator -S $SEED --iso-pep -Y \
    -w "$DATA_OUTPUT_DIR/isotonic/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/isotonic/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/isotonic/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/isotonic/0/peptide.decoy.txt" \
    /data/data/make-pin.pin 2>&1 | tee "$ISOTONIC_NORMAL_DIR/percolator.log"

echo "[INFO]END: $(date)"