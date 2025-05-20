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
INTERM_DIR="$WORK_DIR/preparation/$3/run$4/interm"
MASTER_DIR="$WORK_DIR/$OUTPUT_DIR/irls"
QPAVA_DIR="$WORK_DIR/$OUTPUT_DIR/pava.rank"
QISPLINE_DIR="$WORK_DIR/$OUTPUT_DIR/ispline.rank"
DPAVA_DIR="$WORK_DIR/$OUTPUT_DIR/pava.score"
DISPLINE_DIR="$WORK_DIR/$OUTPUT_DIR/ispline.score"
DATA_OUTPUT_DIR="/data/$OUTPUT_DIR"
DATA_INTERM_DIR="/data/preparation/$3/run$4/interm"

cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }
# Create necessary directories
mkdir -p "$2" "$2/$3" "$OUTPUT_DIR" "$QPAVA_DIR" "$QISPLINE_DIR" "$DPAVA_DIR" "$DISPLINE_DIR"

cp -r "$WORK_DIR/preparation/$3/run$4/irls" "$WORK_DIR/$OUTPUT_DIR"

for KNOCKING_OUT_NUM in 0 "${KNOCKING_OUT_NUM_LIST[@]}"; do
    CURRENT_MASTER_DIR="$MASTER_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QPAVA_DIR="$QPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QISPLINE_DIR="$QISPLINE_DIR/$KNOCKING_OUT_NUM"
    CURRENT_DPAVA_DIR="$DPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_DISPLINE_DIR="$DISPLINE_DIR/$KNOCKING_OUT_NUM"
    mkdir -p "$CURRENT_MASTER_DIR" "$CURRENT_QPAVA_DIR" "$CURRENT_QISPLINE_DIR" "$CURRENT_DPAVA_DIR" "$CURRENT_DISPLINE_DIR"
done

apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
    /usr/bin/percolator -S "$SEED" --pava-pep -Y \
    -w "$DATA_OUTPUT_DIR/pava.rank/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/pava.rank/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/pava.rank/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/pava.rank/0/peptide.decoy.txt" \
    "$DATA_INTERM_DIR/make-pin.pin" 2>&1 | tee "$QPAVA_DIR/0/percolator.log"

apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
    /usr/bin/percolator -S "$SEED" -Y \
    -w "$DATA_OUTPUT_DIR/ispline.rank/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/ispline.rank/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/ispline.rank/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/ispline.rank/0/peptide.decoy.txt" \
    "$DATA_INTERM_DIR/make-pin.pin" 2>&1 | tee "$QISPLINE_DIR/0/percolator.log"

apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
    /usr/bin/percolator -S "$SEED" --pava-pep --ip-pep -Y \
    -w "$DATA_OUTPUT_DIR/pava.score/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/pava.score/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/pava.score/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/pava.score/0/peptide.decoy.txt" \
    "$DATA_INTERM_DIR/make-pin.pin" 2>&1 | tee "$DPAVA_DIR/0/percolator.log"

apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
    /usr/bin/percolator -S "$SEED" --ip-pep -Y \
    -w "$DATA_OUTPUT_DIR/ispline.score/0/weights.pin" \
    -J "$DATA_OUTPUT_DIR/ispline.score/0/features.pin" \
    -r "$DATA_OUTPUT_DIR/ispline.score/0/peptide.target.txt" \
    -B "$DATA_OUTPUT_DIR/ispline.score/0/peptide.decoy.txt" \
    "$DATA_INTERM_DIR/make-pin.pin" 2>&1 | tee "$DISPLINE_DIR/0/percolator.log"

# Loop through the knocking out numbers for running Percolator
for KNOCKING_OUT_NUM in "${KNOCKING_OUT_NUM_LIST[@]}"; do
    CURRENT_MASTER_DIR="$MASTER_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QPAVA_DIR="$QPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_QISPLINE_DIR="$QISPLINE_DIR/$KNOCKING_OUT_NUM"
    CURRENT_DPAVA_DIR="$DPAVA_DIR/$KNOCKING_OUT_NUM"
    CURRENT_DISPLINE_DIR="$DISPLINE_DIR/$KNOCKING_OUT_NUM"

    # Run Percolator with spline regression
    echo "[INFO]Running non-isotonic Percolator using spline regression..."
    apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
        /usr/bin/percolator -S "$SEED" --irls-pep -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/irls/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/irls/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_INTERM_DIR/features.$KNOCKING_OUT_NUM.pin" 2>&1 | tee "$CURRENT_MASTER_DIR/percolator.log"

    # Run Percolator using q-value-based PAVA regression
    echo "[INFO]Running Percolator using q-value-based PAVA regression..."
    apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
        /usr/bin/percolator -S "$SEED" --pava-pep -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/pava.rank/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/pava.rank/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_INTERM_DIR/features.$KNOCKING_OUT_NUM.pin" \
        2>&1 | tee "$CURRENT_QPAVA_DIR/percolator.log"

    # Run Percolator using q-value-based I-Spline regression
    echo "[INFO]Running Percolator using q-value-based I-Spline regression..."
    apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
        /usr/bin/percolator -S "$SEED" -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/ispline.rank/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/ispline.rank/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_INTERM_DIR/features.$KNOCKING_OUT_NUM.pin" \
        2>&1 | tee "$CURRENT_QISPLINE_DIR/percolator.log"
    
    # Run Percolator using decoy probability-based PAVA regression
    echo "[INFO]Running Percolator using decoy probability-based PAVA regression..."
    apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
        /usr/bin/percolator -S "$SEED" --pava-pep --ip-pep -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/pava.score/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/pava.score/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_INTERM_DIR/features.$KNOCKING_OUT_NUM.pin" \
        2>&1 | tee "$CURRENT_DPAVA_DIR/percolator.log"
    
    # Run Percolator using decoy probability-based I-Spline regression
    echo "[INFO]Running Percolator using decoy probability-based I-Spline regression..."
    apptainer run -B "$WORK_DIR:/data" "${WORK_DIR}/percolator.sif" \
        /usr/bin/percolator -S "$SEED" --ip-pep -i 0 -Y \
        -r "$DATA_OUTPUT_DIR/ispline.score/$KNOCKING_OUT_NUM/peptide.target.txt" \
        -B "$DATA_OUTPUT_DIR/ispline.score/$KNOCKING_OUT_NUM/peptide.decoy.txt" \
        "$DATA_INTERM_DIR/features.$KNOCKING_OUT_NUM.pin" \
        2>&1 | tee "$CURRENT_DISPLINE_DIR/percolator.log"
    
done

echo "[INFO]Done."
echo "[INFO]END: $(date)"
