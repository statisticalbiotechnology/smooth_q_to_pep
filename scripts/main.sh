#!/bin/bash
set -e

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <KNOCKING_OUT_NUM_LIST> [SEED] [OUT_DIR]"
    exit 1
fi

echo "[INFO]START: $(date)"
KNOCKING_OUT_NUM_LIST=(${1//[ ,]/ }) # Convert comma or space-separated list to array
SEED=1
OUTPUT_DIR="output"
if [ ! -z "$2" ]; then
    SEED=$2
fi
if [ ! -z "$3" ]; then
    OUTPUT_DIR=$3
fi

WORK_DIR=$(realpath "$PWD/..")

# Define base directories
SCRIPTS_DIR="$WORK_DIR/scripts"
MASTER_DIR="$WORK_DIR/$OUTPUT_DIR/splinePEP"
ISOTONIC_DIR="$WORK_DIR/$OUTPUT_DIR/isoPEP"
IP_ISOTONIC_DIR="$WORK_DIR/$OUTPUT_DIR/ipPEP"
TDC_DIR="$WORK_DIR/$OUTPUT_DIR/tdcPEP"
MASTER_NORMAL_DIR="$MASTER_DIR/0"
REVISED_FEATURES_DIR="$MASTER_NORMAL_DIR/revised_features_pins"
DATA_OUTPUT_DIR="/data/$OUTPUT_DIR"

cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }
# Create necessary directories
mkdir -p "$OUTPUT_DIR" "$MASTER_DIR" "$ISOTONIC_DIR" "$IP_ISOTONIC_DIR" "$TDC_DIR" "$MASTER_NORMAL_DIR" "$REVISED_FEATURES_DIR"

# Run Percolator with spline regression
if [ ! -f "$REVISED_FEATURES_DIR/inserted.features.pin" ]; then
    echo "[INFO]Running non-isotonic Percolator with spline regression..."
    podman run --rm -v "$WORK_DIR:/data" percolator:master \
        /usr/bin/percolator -S $SEED --spline-pep -Y \
        -w "$DATA_OUTPUT_DIR/splinePEP/0/weights.pin" \
        -J "$DATA_OUTPUT_DIR/splinePEP/0/features.pin" \
        -r "$DATA_OUTPUT_DIR/splinePEP/0/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/splinePEP/0/peptide.decoy.txt" \
        /data/data/make-pin.pin 2>&1 | tee "$MASTER_NORMAL_DIR/percolator.log"

    # Insert default direction
    python "$SCRIPTS_DIR/insert_default_direction.py" \
        --weights_pin "$MASTER_NORMAL_DIR/weights.pin" \
        --features_pin "$WORK_DIR/data/make-pin.pin" \
        --output_pin "$REVISED_FEATURES_DIR/inserted.features.pin"
else
    echo "[INFO] Skipping running Percolator with spline regression and inserting default directions as $REVISED_FEATURES_DIR/inserted.features.pin already exists..."
fi
# Loop through the knocking out numbers
for KNOCKING_OUT_NUM in "${KNOCKING_OUT_NUM_LIST[@]}"; do
    echo "[INFO]Changing the label of $KNOCKING_OUT_NUM top-scoring target PSMs..."

    # Define directories for this iteration
    CURRENT_MASTER_DIR="$MASTER_DIR/$KNOCKING_OUT_NUM"
    CURRENT_ISOTONIC_DIR="$ISOTONIC_DIR/$KNOCKING_OUT_NUM"
    CURRENT_IP_ISOTONIC_DIR="$IP_ISOTONIC_DIR/$KNOCKING_OUT_NUM"
    CURRENT_TDC_DIR="$TDC_DIR/$KNOCKING_OUT_NUM"
    mkdir -p "$CURRENT_MASTER_DIR" "$CURRENT_ISOTONIC_DIR" "$CURRENT_IP_ISOTONIC_DIR" "$CURRENT_TDC_DIR"

    # Revise labels
    python "$SCRIPTS_DIR/revise_labels.py" \
        --target_peptide_file "$MASTER_NORMAL_DIR/peptide.target.txt" \
        --features_file "$REVISED_FEATURES_DIR/inserted.features.pin" \
        --revised_features_file "$REVISED_FEATURES_DIR/features.$KNOCKING_OUT_NUM.pin" \
        --N "$KNOCKING_OUT_NUM"

    # Run standard isotonic Percolator with default tdc_to_pep
    echo "[INFO]Running isotonic Percolator using TDC..."
    podman run --rm -v "$WORK_DIR:/data" percolator:master \
        /usr/bin/percolator -S $SEED -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/tdcPEP/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/tdcPEP/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_OUTPUT_DIR/splinePEP/0/revised_features_pins/features.$KNOCKING_OUT_NUM.pin" 2>&1 | tee "$CURRENT_ISOTONIC_DIR/percolator.log"

    # Run isotonic Percolator with --pep-from-q
    echo "[INFO]Running isotonic Percolator using isotonic regression..."
    podman run --rm -v "$WORK_DIR:/data" percolator:master \
        /usr/bin/percolator -S $SEED --pep-from-q -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/isoPEP/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/isoPEP/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_OUTPUT_DIR/splinePEP/0/revised_features_pins/features.$KNOCKING_OUT_NUM.pin" 2>&1 | tee "$CURRENT_ISOTONIC_DIR/percolator.log"

    # Run isotonic Percolator with --ip-pep
    echo "[INFO]Running isotonic Percolator using isotonic regression with interpolation..."
    podman run --rm -v "$WORK_DIR:/data" percolator:master \
        /usr/bin/percolator -S $SEED --ip-pep -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/ipPEP/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/ipPEP/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_OUTPUT_DIR/splinePEP/0/revised_features_pins/features.$KNOCKING_OUT_NUM.pin" 2>&1 | tee "$CURRENT_IP_ISOTONIC_DIR/percolator.log"

    # Run Percolator with spline regression
    echo "[INFO]Running non-isotonic Percolator using spline regression..."
    podman run --rm -v "$WORK_DIR:/data" percolator:master \
        /usr/bin/percolator -S $SEED --spline-pep -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/splinePEP/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/splinePEP/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_OUTPUT_DIR/splinePEP/0/revised_features_pins/features.$KNOCKING_OUT_NUM.pin" 2>&1 | tee "$CURRENT_MASTER_DIR/percolator.log"
done

echo "[INFO]Done."
echo "[INFO]END: $(date)"
