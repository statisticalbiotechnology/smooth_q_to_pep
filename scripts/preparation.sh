#!/bin/bash
set -e

# On Tetralith: /proj/proteoforma_nsc/smooth_q_to_pep
# apptainer build percolator.sif docker://ghcr.io/statisticalbiotechnology/percolator:isotonic
# apptainer build pyisopep.sif docker://ghcr.io/statisticalbiotechnology/pyisotonicpep:main
# apptainer run --pwd /app pyisopep.sif --help

# Check four arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <KNOCKING_OUT_NUM_LIST> <OUT_DIR> <DATASET> <SEED>"
    exit 1
fi

echo "[INFO]START: $(date)"
KNOCKING_OUT_NUM_LIST=(${1//[ ,]/ }) # Convert comma or space-separated list to array
SEED=$4
OUTPUT_DIR="$2/$3/run$4"

# WORK_DIR=$(realpath "$PWD/..")  # /proj/proteoforma_nsc/smooth_q_to_pep
WORK_DIR="/proj/proteoforma_nsc/smooth_q_to_pep"

# Define base directories
SCRIPTS_DIR="$WORK_DIR/scripts"
INTERM_DIR="$WORK_DIR/$OUTPUT_DIR/interm"
MASTER_DIR="$WORK_DIR/$OUTPUT_DIR/irls"
DATA_OUTPUT_DIR="/data/$OUTPUT_DIR"
DATA_INTERM_DIR="/data/$OUTPUT_DIR/interm"

cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }
# Create necessary directories
mkdir -p "$2" "$2/$3" "$OUTPUT_DIR" "$INTERM_DIR" "$MASTER_DIR" "$MASTER_DIR/0"

# Run crux
if [ ! -f "$INTERM_DIR/make-pin.pin" ]; then
    echo "[INFO]Running Crux..."
    bash "$WORK_DIR/datasets/$3/run_crux.sh" "$2" "$SEED"
else
    echo "[INFO] Skipping running Crux as $INTERM_DIR/make-pin.pin already exists..."
fi

# Run Percolator with spline regression
if [ ! -f "$INTERM_DIR/features.0.pin" ]; then
    echo "[INFO]Running non-isotonic Percolator with spline regression..."
    apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
        /usr/bin/percolator -S "$SEED" --irls-pep -Y \
        -w "$DATA_OUTPUT_DIR/irls/0/weights.pin" \
        -J "$DATA_OUTPUT_DIR/irls/0/features.pin" \
        -r "$DATA_OUTPUT_DIR/irls/0/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/irls/0/peptide.decoy.txt" \
        "$DATA_INTERM_DIR/make-pin.pin" 2>&1 | tee "$MASTER_DIR/0/percolator.log"

    cp "$MASTER_DIR/0/peptide.target.txt" "$INTERM_DIR/peptide.target.0.txt"

    # Insert default direction
    python "$SCRIPTS_DIR/insert_default_direction.py" \
        --weights_pin "$MASTER_DIR/0/weights.pin" \
        --features_pin "$INTERM_DIR/make-pin.pin" \
        --output_pin "$INTERM_DIR/features.0.pin"
else
    echo "[INFO] Skipping running Percolator with spline regression and inserting default directions as $INTERM_DIR/features.0.pin already exists..."
fi

# Loop through the knocking out numbers for revising labels and running Percolator
for KNOCKING_OUT_NUM in "${KNOCKING_OUT_NUM_LIST[@]}"; do
    echo "[INFO]Changing the label of $KNOCKING_OUT_NUM top-scoring target PSMs..."
    # Revise labels for the inserted.features.pin file
    python "$SCRIPTS_DIR/revise_labels.py" \
        --target_peptide_file "$INTERM_DIR/peptide.target.0.txt" \
        --features_file "$INTERM_DIR/features.0.pin" \
        --revised_features_file "$INTERM_DIR/features.$KNOCKING_OUT_NUM.pin" \
        --N "$KNOCKING_OUT_NUM"
done

echo "[INFO]Done."
echo "[INFO]END: $(date)"
