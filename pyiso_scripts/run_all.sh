#!/bin/bash
set -e

# Define the list of PXD codes
PXD_CODES=(
#     "PXD003868"
#     "PXD004325"
    "PXD004424"
    "PXD004467"
    "PXD004536"
    "PXD004565"
    "PXD004947"
    "PXD004948"
    "PXD005025"
    "PXD013274"
)

# Loop over each code and run repeat_main.sh
for code in "${PXD_CODES[@]}"; do
    echo "[INFO] Running repeat_main.sh for pyiso_test/$code"
    ./repeat_main.sh 10 "1,2,3,4" "pyiso_test/$code"
done

echo "[INFO] All jobs finished."