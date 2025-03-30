#!/bin/bash
set -e

# On Tetralith: /proj/proteoforma_nsc/smooth_q_to_pep/pyiso_scripts
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
SCRIPTS_DIR="$WORK_DIR/pyiso_scripts"
INTERM_DIR="$WORK_DIR/$OUTPUT_DIR/interm"
MASTER_DIR="$WORK_DIR/$OUTPUT_DIR/percolator"
QPAVA_DIR="$WORK_DIR/$OUTPUT_DIR/qPAVA"
QIPPAVA_DIR="$WORK_DIR/$OUTPUT_DIR/qipPAVA"
QISPLINE_DIR="$WORK_DIR/$OUTPUT_DIR/qISpline"
DPAVA_DIR="$WORK_DIR/$OUTPUT_DIR/dPAVA"
DISPLINE_DIR="$WORK_DIR/$OUTPUT_DIR/dISpline"
DATA_OUTPUT_DIR="/data/$OUTPUT_DIR"
DATA_INTERM_DIR="/data/$OUTPUT_DIR/interm"

cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }
# Create necessary directories
mkdir -p "$2" "$2/$3" "$OUTPUT_DIR" "$INTERM_DIR" "$MASTER_DIR" "$MASTER_DIR/0" "$QPAVA_DIR" "$QIPPAVA_DIR" "$QISPLINE_DIR" "$DPAVA_DIR" "$DISPLINE_DIR"

# Run crux
if [ ! -f "$INTERM_DIR/make-pin.pin" ]; then
    echo "[INFO]Running Crux..."
    bash "$WORK_DIR/datasets/$3/run_crux.sh" "$2" "$SEED"
else
    echo "[INFO] Skipping running Crux as $INTERM_DIR/make-pin.pin already exists..."
fi

# Run Percolator with spline regression
if [ ! -f "$INTERM_DIR/inserted.features.pin" ]; then
    echo "[INFO]Running non-isotonic Percolator with spline regression..."
    apptainer run -B "$WORK_DIR:/data" "${SCRIPTS_DIR}/percolator.sif" \
        /usr/bin/percolator -S "$SEED" --spline-pep -Y \
        -w "$DATA_OUTPUT_DIR/percolator/0/weights.pin" \
        -J "$DATA_OUTPUT_DIR/percolator/0/features.pin" \
        -r "$DATA_OUTPUT_DIR/percolator/0/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/percolator/0/peptide.decoy.txt" \
        "$DATA_INTERM_DIR/make-pin.pin" 2>&1 | tee "$MASTER_DIR/0/percolator.log"
    
    cp "$MASTER_DIR/0/peptide.target.txt" "$INTERM_DIR/peptide.target.0.txt"
    cp "$MASTER_DIR/0/peptide.decoy.txt" "$INTERM_DIR/peptide.decoy.0.txt"

    # Insert default direction
    python "$SCRIPTS_DIR/insert_default_direction.py" \
        --weights_pin "$MASTER_DIR/0/weights.pin" \
        --features_pin "$INTERM_DIR/make-pin.pin" \
        --output_pin "$INTERM_DIR/inserted.features.pin"
else
    echo "[INFO] Skipping running Percolator with spline regression and inserting default directions as $INTERM_DIR/inserted.features.pin already exists..."
fi

for KNOCKING_OUT_NUM in 0 "${KNOCKING_OUT_NUM_LIST[@]}"; do
    CURRENT_MASTER_DIR="$MASTER_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QPAVA_DIR="$QPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QIPPAVA_DIR="$QIPPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QISPLINE_DIR="$QISPLINE_DIR/$KNOCKING_OUT_NUM"
    CURRENT_DPAVA_DIR="$DPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_DISPLINE_DIR="$DISPLINE_DIR/$KNOCKING_OUT_NUM"
    mkdir -p "$CURRENT_MASTER_DIR" "$CURRENT_QPAVA_DIR" "$CURRENT_QIPPAVA_DIR" "$CURRENT_QISPLINE_DIR" "$CURRENT_DPAVA_DIR" "$CURRENT_DISPLINE_DIR"
done

# Loop through the knocking out numbers for revising labels and running Percolator
for KNOCKING_OUT_NUM in "${KNOCKING_OUT_NUM_LIST[@]}"; do
    echo "[INFO]Changing the label of $KNOCKING_OUT_NUM top-scoring target PSMs..."
    CURRENT_MASTER_DIR="$MASTER_DIR/$KNOCKING_OUT_NUM"
    # Revise labels for the inserted.features.pin file
    python "$SCRIPTS_DIR/revise_labels.py" \
        --target_peptide_file "$INTERM_DIR/peptide.target.0.txt" \
        --features_file "$INTERM_DIR/inserted.features.pin" \
        --revised_features_file "$INTERM_DIR/features.$KNOCKING_OUT_NUM.pin" \
        --N "$KNOCKING_OUT_NUM"
    # Revise the target and decoy files
    python "$SCRIPTS_DIR/revise_files.py" \
        --target "$INTERM_DIR/peptide.target.0.txt" \
        --decoy "$INTERM_DIR/peptide.decoy.0.txt" \
        --revised_target "$INTERM_DIR/peptide.target.$KNOCKING_OUT_NUM.txt" \
        --revised_decoy "$INTERM_DIR/peptide.decoy.$KNOCKING_OUT_NUM.txt" \
        --N "$KNOCKING_OUT_NUM"
    
    # Run Percolator with spline regression
    echo "[INFO]Running non-isotonic Percolator using spline regression..."
    apptainer run -B "$WORK_DIR:/data" "${SCRIPTS_DIR}/percolator.sif" \
        /usr/bin/percolator -S "$SEED" --spline-pep -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/percolator/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/percolator/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_INTERM_DIR/features.$KNOCKING_OUT_NUM.pin" 2>&1 | tee "$CURRENT_MASTER_DIR/percolator.log"
done

# Loop through 0 and the knocking out numbers for running pyIsotonicPEP
for KNOCKING_OUT_NUM in 0 "${KNOCKING_OUT_NUM_LIST[@]}"; do
    CURRENT_QPAVA_DIR="$QPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QIPPAVA_DIR="$QIPPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QISPLINE_DIR="$QISPLINE_DIR/$KNOCKING_OUT_NUM"
    CURRENT_DPAVA_DIR="$DPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_DISPLINE_DIR="$DISPLINE_DIR/$KNOCKING_OUT_NUM"
    # Run pyIsotonicPEP using q-value-based PAVA regression
    echo "[INFO]Running pyIsotonicPEP using q-value-based PAVA regression..."
    apptainer run --pwd /app -B "$WORK_DIR:/data" "${SCRIPTS_DIR}/pyisopep.sif" \
        q2pep --input "$DATA_INTERM_DIR/peptide.target.$KNOCKING_OUT_NUM.txt" \
        --regression-algo "PAVA" \
        --output "$DATA_OUTPUT_DIR/qPAVA/$KNOCKING_OUT_NUM/peptide.target.txt" \
        2>&1 | tee "$CURRENT_QPAVA_DIR/pyisopep.log"

    # Run pyIsotonicPEP using q-value-based PAVA interpolation
    echo "[INFO]Running pyIsotonicPEP using q-value-based PAVA interpolation..."
    apptainer run --pwd /app -B "$WORK_DIR:/data" "${SCRIPTS_DIR}/pyisopep.sif" \
        q2pep --input "$DATA_INTERM_DIR/peptide.target.$KNOCKING_OUT_NUM.txt" \
        --regression-algo "PAVA" --ip \
        --output "$DATA_OUTPUT_DIR/qipPAVA/$KNOCKING_OUT_NUM/peptide.target.txt" \
        2>&1 | tee "$CURRENT_QIPPAVA_DIR/pyisopep.log"

    # Run pyIsotonicPEP using q-value-based I-Spline regression
    echo "[INFO]Running pyIsotonicPEP using q-value-based I-Spline regression..."
    apptainer run --pwd /app -B "$WORK_DIR:/data" "${SCRIPTS_DIR}/pyisopep.sif" \
        q2pep --input "$DATA_INTERM_DIR/peptide.target.$KNOCKING_OUT_NUM.txt" \
        --output "$DATA_OUTPUT_DIR/qISpline/$KNOCKING_OUT_NUM/peptide.target.txt" \
        2>&1 | tee "$CURRENT_QISPLINE_DIR/pyisopep.log"
    
    # Run pyIsotonicPEP using decoy probability-based I-Spline regression
    echo "[INFO]Running pyIsotonicPEP using decoy probability-based PAVA regression..."
    apptainer run --pwd /app -B "$WORK_DIR:/data" "${SCRIPTS_DIR}/pyisopep.sif" \
        d2pep --target-file "$DATA_INTERM_DIR/peptide.target.$KNOCKING_OUT_NUM.txt" \
        --decoy-file "$DATA_INTERM_DIR/peptide.decoy.$KNOCKING_OUT_NUM.txt" \
        --regression-algo "PAVA" \
        --output "$DATA_OUTPUT_DIR/dPAVA/$KNOCKING_OUT_NUM/peptide.target.txt" \
        2>&1 | tee "$CURRENT_DPAVA_DIR/pyisopep.log"
    
    # Run pyIsotonicPEP using decoy probability-based I-Spline regression
    echo "[INFO]Running pyIsotonicPEP using decoy probability-based I-Spline regression..."
    apptainer run --pwd /app -B "$WORK_DIR:/data" "${SCRIPTS_DIR}/pyisopep.sif" \
        d2pep --target-file "$DATA_INTERM_DIR/peptide.target.$KNOCKING_OUT_NUM.txt" \
        --decoy-file "$DATA_INTERM_DIR/peptide.decoy.$KNOCKING_OUT_NUM.txt" \
        --output "$DATA_OUTPUT_DIR/dISpline/$KNOCKING_OUT_NUM/peptide.target.txt" \
        2>&1 | tee "$CURRENT_DISPLINE_DIR/pyisopep.log"
    
done

echo "[INFO]Done."
echo "[INFO]END: $(date)"
