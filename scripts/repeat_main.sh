#!/bin/bash
set -e

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <NUM_RUNS> <KNOCKING_OUT_NUM_LIST> <OUT_DIR>"
    exit 1
fi

NUM_RUNS=$1
KNOCKING_OUT_NUM_LIST=$2
OUT_DIR=$3
MAIN_SCRIPT="./main.sh"
ISOPERC_SCRIPT="./run_normal_isoperc.sh"

if [ ! -f "$MAIN_SCRIPT" ]; then
    echo "[ERROR] $MAIN_SCRIPT not found!"
    exit 1
fi

if [ ! -f "$ISOPERC_SCRIPT" ]; then
    echo "[ERROR] $ISOPERC_SCRIPT not found!"
    exit 1
fi

for ((i=1; i<=NUM_RUNS; i++)); do
    echo "[INFO] Running run $i with SEED=$i and OUTPUT_DIR=$OUT_DIR..."
    RUN_DIR="$OUT_DIR/run$i"
    SEED=$i
    bash $MAIN_SCRIPT "$KNOCKING_OUT_NUM_LIST" $SEED $RUN_DIR
    bash $ISOPERC_SCRIPT $SEED $RUN_DIR

done

echo "[INFO] All iterations completed successfully."