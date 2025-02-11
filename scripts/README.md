# README: Running Percolator with and without Isotonic Regression

This document provides instructions on how to run Percolator with and without isotonic regression.

---

## Preparation

1. **Generate PIN TSV File (if needed):**
   If a PIN TSV file has not been generated yet, use the following commands to create it using the Crux toolkit:
   ```bash
   cd scripts
   ./run_crux.sh
   ```

---

## Running Percolator

2. **Run Percolator with and without Isotonic Regression:**
   Use the following script to run Percolator with the specified knocking-out numbers list:
   ```bash
   cd scripts
   ./main.sh <KNOCKING_OUT_NUM_LIST>
   ```
   Example:
   ```bash
   ./main.sh "1,5,10,100" 2>&1 | tee logs/main_$(date '+%Y%m%d_%H%M%S').log
   ```

---

## Notes

- Ensure that the `scripts` directory contains all the necessary scripts, including `run_crux.sh` and `main.sh`.
- Create a `logs` directory beforehand to store log files:
  ```bash
  mkdir -p logs
  ```
- The `<KNOCKING_OUT_NUM_LIST>` can be provided as a space-separated or comma-separated list of numbers (e.g., `1,5,10,100` or `1 5 10 100`).

---

## Example Workflow

```bash
# Step 1: Generate PIN TSV File
cd scripts
./run_crux.sh

# Step 2: Run Percolator and Save Logs
./main.sh "1,5,10,100" 2>&1 | tee logs/main_$(date '+%Y%m%d_%H%M%S').log

# [Optional] Step 3: Run isotonic Percolator for comparison (if needed).
# This step is optional and not mandatory for the primary workflow. Use it if you want to compare normal Percolator runs with and without isotonic regression.
./run_normal_isoperc.sh
```

### Repeat

```bash
./repeat_main.sh <NUM_RUNS> <KNOCKING_OUT_NUM_LIST> <OUT_DIR>
./repeat_main.sh 5 "1,2,3,4" "output_rep"
```

---
