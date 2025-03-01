# Isotonic PEP Estimation

## Overview
[**pyIsotonicPEP**](https://github.com/statisticalbiotechnology/smooth_q_to_pep/pkgs/container/pyisotonicpep/365477588?tag=main) provides a unified interface for estimating Posterior Error Probabilities (PEPs) using isotonic regression for identifications in shotgun proteomics. It supports two methods:
- **q2pep:** Estimate non-decreasing PEPs from q-values.
- **obs2pep:** Estimate non-decreasing PEPs from a stream of target and decoy observations derived from target-decoy competition (TDC) method.

The package consists of two main Python files:
- [**`IsotonicPEP.py`**](https://github.com/statisticalbiotechnology/smooth_q_to_pep/blob/main/pyIsoPEP/IsotonicPEP.py): Implements isotonic regression and methods for PEP estimation.
- [**`main.py`**](https://github.com/statisticalbiotechnology/smooth_q_to_pep/blob/main/pyIsoPEP/main.py): Provides a command-line interface (CLI) to run the isotonic regressor with various options.

## Features
- **Isotonic Regression:** Enforces a non-decreasing constraint on probability estimates.
- **q2pep Method:** Converts q-values to PEP values.
- **obs2pep Method:** Processes target-decoy observations to compute PEP values.
- **Optional PEP-based q-value Estimation:** Calculates q-values from the estimated PEPs.
- **Flexible Input Formats:** For the obs2pep Method, supports both concatenated input files and separate target/decoy files.

## Usage
### Pull the docker image and run the container to display the help message
```bash
$ podman pull ghcr.io/statisticalbiotechnology/pyisotonicpep:main
$ podman run --rm -it pyisotonicpep:main -h
```
### Command-Line Options
#### General Options
```less
usage: main.py [-h] [--no-calc-q] [--verbose] {q2pep,obs2pep} ...
```
*  **`{q2pep,obs2pep}`**    Select the PEP estimation method.
    - **`q2pep`** Estimate PEP from q-values (q2pep method).
    - **`obs2pep`**   Estimate PEP from target-decoy observations(obs2pep method).
*  **`-h, --help`** Show this help message and exit.
*  **`--no-calc-q`**    Do not estimate q-values from calculated PEPs (default: calculate q-values).
*  **`--verbose`**  Print parameter information.
#### q2pep Mode Options
```less
usage: main.py q2pep [-h] --input INPUT [--qcol QCOL] [--smooth] [--pava-method {basic,ip}] [--center-method {mean,median}] --output OUTPUT
```
*  **`--input INPUT`**  Path to the TSV file containing q-values.
*  **`--qcol QCOL`**    Column name for q-values in the input file (default: 'q-value').
*  **`--smooth`**   Apply block-merge pre-processing (default: False).
    *  **`--pava-method {basic,ip}`**   Choose PAVA method to use: 'basic' for  PAVA regression or 'ip' for PAVA interpolation (default: basic).
*  **`--center-method {mean,median}`**  ChoosePAVA center method to use in PAVA interpolation (default: mean).
*  **`--output OUTPUT`**    Output file path (or output directory; if directory, a default name 'outputPEP.target.txt' will be used).
#### obs2pep  Mode Options
```less
usage: main.py obs2pep [-h] (--cat-file CAT_FILE | --target-file TARGET_FILE) [--decoy-file DECOY_FILE] [--score-col SCORE_COL] [--type-col TYPE_COL] [--target-label TARGET_LABEL] [--decoy-label DECOY_LABEL] -output OUTPUT
```
*  **`--cat-file CAT_FILE`**   Path to a concatenated TSV file containing score and label columns.
*  **`--target-file TARGET_FILE`**  Path to the TSV file containing target scores (used in separate input mode).
*  **`--decoy-file DECOY_FILE`**    Path to the TSV file containing decoy scores (required in separate input mode).
*  **`--score-col SCORE_COL`**  Column name for score (default: 'score').
*  **`--type-col TYPE_COL`**    Column name for target/decoy (default: 'label').
*  **`--target-label TARGET_LABEL`**    Target label in concatenated input file (default: 'target').
*  **`--decoy-label DECOY_LABEL`**  Decoy label in concatenated input file (default: 'decoy').
*  **`--output OUTPUT`**    Output file path or directory for target results. (Note: only target PEPs are saved.)

### Running examples with the provided files
#### 1. q2pep Mode
```bash
$ podman run --rm -it -v .:/data pyisotonicpep:main q2pep --input example/peptide.target.txt --output /data
```
#### 2. obs2pep Mode
**2.1 using a concatenated target and decoy input**
```bash
$ podman run --rm -it -v .:/data pyisotonicpep:main obs2pep --cat-file example/peptide.cat.txt --score-col score --type-col type --target-label 0 --decoy-label 1 --output /data
```
**2.2 using separate target and decoy inputs**
```bash
$ podman run --rm -it -v .:/data pyisotonicpep:main obs2pep --target-file example/peptide.target.txt --decoy-file example/peptide.decoy.txt --score-col score --output /data
```