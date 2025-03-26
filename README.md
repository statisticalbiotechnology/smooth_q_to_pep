# Isotonic PEP Estimation

## Overview
[**pyIsotonicPEP**](https://github.com/statisticalbiotechnology/smooth_q_to_pep/pkgs/container/pyisotonicpep/382208731?tag=main) provides a unified interface for estimating Posterior Error Probabilities (PEPs) using isotonic regression for identifications in shotgun proteomics. It supports two methods:
- **q2pep:** Estimate non-decreasing PEPs from q-values.
- **d2pep:** Estimate non-decreasing PEPs from a stream of target and decoy observations derived from target-decoy competition (TDC) method.

The package consists of two main Python files:
- [**`IsotonicPEP.py`**](https://github.com/statisticalbiotechnology/smooth_q_to_pep/blob/main/pyIsoPEP/IsotonicPEP.py): Implements isotonic regression for PEP estimation.
- [**`main.py`**](https://github.com/statisticalbiotechnology/smooth_q_to_pep/blob/main/pyIsoPEP/main.py): Provides a command-line interface (CLI) to run the isotonic regressor with various options.

## Features
- **Isotonic Regression:** Enforces a non-decreasing constraint on probability estimates.
- **q2pep Method:** Converts q-values to PEP values.
- **d2pep Method:** Processes target-decoy observations to compute PEP values.
- **Optional PEP-based q-value Estimation:** Calculates q-values from the estimated PEPs.
- **Flexible Input Formats:** In `d2pep` mode, supports both concatenated input files and separate target/decoy files.
- **Configurable Regression Options:** Choose between PAVA and I-Spline regression, with options for interpolation and block center calculation.


## Usage
### Pull the Docker image from GitHub Container Registry and display the help message
```bash
$ podman pull ghcr.io/statisticalbiotechnology/pyisotonicpep:main
$ podman run --rm -it pyisotonicpep:main -h
```
### Command-Line Options
#### General Options
```less
usage: main.py [-h] [--no-calc-q] [--verbose] {q2pep,d2pep} ...
```
*  **`{q2pep,d2pep}`**    Select the PEP estimation method.
    - **`q2pep`** Estimate PEPs from q-values.
    - **`d2pep`**   Estimate PEPs from target-decoy observations.
*  **`-h, --help`** Show this help message and exit.
*  **`--no-calc-q`**    Do not estimate q-values from calculated PEPs (by default, q-values are calculated).
*  **`--verbose`**  Print detailed parameter information.
#### q2pep Mode Options
```less
usage: main.py q2pep [-h] --input INPUT [--qcol QCOL] [--regression-algo {PAVA,ispline}] [--max-iter MAX_ITER] [--ip] [--ip-algo {ispline,pchip}] [--center-method {mean,median}] --output OUTPUT
```
*  **`--input INPUT`**  Path to the TSV file containing q-values.
*  **`--qcol QCOL`**    Column name for q-values in the input file (default: 'q-value').
*  **`--regression-algo {PAVA,ispline}`**  Regression algorithm to use (default: 'ispline').
*  **`--max-iter MAX_ITER`**  Maximum iterations for I-Spline iterative solvers (default: 5000).
*  **`--ip`**  Apply monotonic interpolation (only used when using the PAVA regression algorithm).
*  **`--ip-algo {ispline,pchip}`**  Interpolation algorithm for PAVA-derived block centers (default: 'ispline').
*  **`--center-method {mean,median}`**  Method for computing block centers in interpolation (default: 'mean').
*  **`--output OUTPUT`**  Output file path or directory. If a directory is provided, the default filename outputPEP.target.qbased.txt will be used.  
#### d2pep  Mode Options
```less
usage: main.py d2pep [-h] (--cat-file CAT_FILE | --target-file TARGET_FILE) [--decoy-file DECOY_FILE] [--score-col SCORE_COL] [--type-col TYPE_COL] [--target-label TARGET_LABEL] [--decoy-label DECOY_LABEL] [--regression-algo {PAVA,ispline}] [--max-iter MAX_ITER] --output OUTPUT
```
*  **`--cat-file CAT_FILE`**   Path to a concatenated TSV file containing score and label columns.
*  **`--target-file TARGET_FILE`**  Path to the TSV file containing target scores (used in separate input mode).
*  **`--decoy-file DECOY_FILE`**    Path to the TSV file containing decoy scores (required in separate input mode).
*  **`--score-col SCORE_COL`**  Column name for score (default: 'score').
*  **`--type-col TYPE_COL`**    Column name for target/decoy (default: 'label').
*  **`--target-label TARGET_LABEL`**    Target label in concatenated input file (default: 'target').
*  **`--decoy-label DECOY_LABEL`**  Decoy label in concatenated input file (default: 'decoy').
*  **`--regression-algo {PAVA,ispline}`**  Regression algorithm to use (default: 'ispline').
*  **`--max-iter MAX_ITER`**  Maximum iterations for I-Spline iterative solvers (default: 5000).
*  **`--output OUTPUT`**    Output file path or directory. If a directory is provided, the default filename outputPEP.target.dbased.txt will be used (note: only target PEPs are saved).

### Running examples with the provided files
#### 1. q2pep Mode
```bash
$ podman run --rm -it -v .:/data pyisotonicpep:main q2pep --input /example/peptide.target.txt --output /data
```
#### 2. d2pep Mode
**2.1 using a concatenated target and decoy input**
```bash
$ podman run --rm -it -v .:/data pyisotonicpep:main d2pep --cat-file /example/peptide.cat.txt --score-col score --type-col type --target-label 0 --decoy-label 1 --output /data
```
**2.2 using separate target and decoy inputs**
```bash
$ podman run --rm -it -v .:/data pyisotonicpep:main d2pep --target-file /example/peptide.target.txt --decoy-file /example/peptide.decoy.txt --score-col score --output /data
```