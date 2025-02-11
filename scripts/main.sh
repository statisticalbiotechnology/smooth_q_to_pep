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
MASTER_DIR="$WORK_DIR/$OUTPUT_DIR/master"
ISOTONIC_DIR="$WORK_DIR/$OUTPUT_DIR/isotonic"
MASTER_NORMAL_DIR="$MASTER_DIR/0"
REVISED_FEATURES_DIR="$MASTER_NORMAL_DIR/revised_features_pins"
DATA_OUTPUT_DIR="/data/$OUTPUT_DIR"

cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }
# Create necessary directories
mkdir -p "$OUTPUT_DIR" "$MASTER_DIR" "$ISOTONIC_DIR" "$MASTER_NORMAL_DIR" "$REVISED_FEATURES_DIR"

# Run Percolator normally
if [ ! -f "$REVISED_FEATURES_DIR/inserted.features.pin" ]; then
    echo "[INFO]Running Percolator normally..."
    podman run --rm -v "$WORK_DIR:/data" percolator:isotonic \
        /usr/bin/percolator -S $SEED -Y \
        -w "$DATA_OUTPUT_DIR/master/0/weights.pin" \
        -J "$DATA_OUTPUT_DIR/master/0/features.pin" \
        -r "$DATA_OUTPUT_DIR/master/0/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/master/0/peptide.decoy.txt" \
        /data/data/make-pin.pin 2>&1 | tee "$MASTER_NORMAL_DIR/percolator.log"

    # Insert default direction
    python "$SCRIPTS_DIR/insert_default_direction.py" \
        --weights_pin "$MASTER_NORMAL_DIR/weights.pin" \
        --features_pin "$WORK_DIR/data/make-pin.pin" \
        --output_pin "$REVISED_FEATURES_DIR/inserted.features.pin"
else
    echo "[INFO] Skipping running Percolator normally and inserting default directions as $REVISED_FEATURES_DIR/inserted.features.pin already exists..."
fi
# Loop through the knocking out numbers
for KNOCKING_OUT_NUM in "${KNOCKING_OUT_NUM_LIST[@]}"; do
    echo "[INFO]Changing the label of $KNOCKING_OUT_NUM top-scoring target PSMs..."

    # Define directories for this iteration
    CURRENT_MASTER_DIR="$MASTER_DIR/$KNOCKING_OUT_NUM"
    CURRENT_ISOTONIC_DIR="$ISOTONIC_DIR/$KNOCKING_OUT_NUM"
    mkdir -p "$CURRENT_MASTER_DIR" "$CURRENT_ISOTONIC_DIR"

    # Revise labels
    python "$SCRIPTS_DIR/revise_labels.py" \
        --target_peptide_file "$MASTER_NORMAL_DIR/peptide.target.txt" \
        --features_file "$REVISED_FEATURES_DIR/inserted.features.pin" \
        --revised_features_file "$REVISED_FEATURES_DIR/features.$KNOCKING_OUT_NUM.pin" \
        --N "$KNOCKING_OUT_NUM"

    # Run isotonic Percolator
    echo "[INFO]Running isotonic Percolator..."
    podman run --rm -v "$WORK_DIR:/data" percolator:isotonic \
        /usr/bin/percolator -S $SEED --iso-pep -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/isotonic/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/isotonic/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_OUTPUT_DIR/master/0/revised_features_pins/features.$KNOCKING_OUT_NUM.pin" 2>&1 | tee "$CURRENT_ISOTONIC_DIR/percolator.log"

    # Run non-isotonic Percolator
    echo "[INFO]Running non-isotonic Percolator..."
    podman run --rm -v "$WORK_DIR:/data" percolator:isotonic \
        /usr/bin/percolator -S $SEED -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/master/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/master/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_OUTPUT_DIR/master/0/revised_features_pins/features.$KNOCKING_OUT_NUM.pin" 2>&1 | tee "$CURRENT_MASTER_DIR/percolator.log"
done

echo "[INFO]Done :)"
echo "[INFO]END: $(date)"
