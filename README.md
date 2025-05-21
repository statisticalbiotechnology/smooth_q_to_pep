# `pyIsoPEP`: Isotonic PEP Estimator

## Overview
`pyIsoPEP` provides both a [**Docker image**](https://github.com/statisticalbiotechnology/smooth_q_to_pep/pkgs/container/pyisotonicpep) and a [**Python API**](https://pypi.org/project/pyIsoPEP/) with a unified interface for estimating smooth, non‑decreasing Posterior Error Probabilities (PEPs) using isotonic regression for target identifications in shotgun proteomics.

Two workflows are available:

| workflow | scheme | starting point |
|----------|----------------|----------------|
| **q2pep** | rank-based | 	a list of target identifications with *q-values* |
| **d2pep** | score-based | target‑decoy competition (TDC) output: a list of target and decoy identifications with *scores* |

Internally, both workflows use isotonic regression - implemented via either the Pool-Adjacent-Violators Algorithm (PAVA) or I‑Splines. Optional post-processing can be applied to derive q-values from the estimated PEPs. In I‑Spline regression, for **d2pep**, scores are adaptively binned to balance decoy counts per bin, with tighter bins at lower scores. Linear weight rescaling is applied to emphasize early (low-PEP) regions.

---

## Installation

### PyPI

```bash
pip install pyIsoPEP
```

### Docker image

```bash
podman pull ghcr.io/statisticalbiotechnology/pyisotonicpep:main
# or
docker pull ghcr.io/statisticalbiotechnology/pyisotonicpep:main
```

---

## Input & output

**Input:** tab‑separated value (TSV) files; column names are configurable.  
**Output:** the original **target** rows plus up to four extra columns:

| column | present when | description |
|--------|--------------|---------|
| `pyIsoPEP FDR` | `--calc-q-from-fdr` | estimated False Discovery Rates (FDRs) |
| `pyIsoPEP q-value from FDR` | `--calc-q-from-fdr` | FDR-derived q-values |
| `pyIsoPEP PEP` | always | estimated PEPs |
| `pyIsoPEP q-value from PEP` | `--calc-q-from-pep` | PEP-derived q-values |

---

## Command-line reference

Top‑level help:

```bash
pyisopep -h
# or
podman run --rm -it pyisotonicpep:main -h
```

Workflow‑specific help:

```bash
pyisopep q2pep -h
pyisopep d2pep -h
# or
podman run --rm -it pyisotonicpep:main q2pep -h
podman run --rm -it pyisotonicpep:main d2pep -h
```

### Common flags (valid for *both* workflows)

| flag | default | description |
|------|---------|---------|
| `--cat-file FILE`, `--target-file FILE` | *required* | Select concatenated vs. separate input |
| `--decoy-file FILE` | — | Required only with separate target & decoy lists |
| `--score-col COL` | `score` | Score column |
| `--label-col COL` | `label` | Column distinguishing targets/decoys in a concatenated file |
| `--target-label STR` | `target` | String marking target rows |
| `--decoy-label STR` | `decoy` | String marking decoy rows |
| `--regression-algo {PAVA,ispline}` | `ispline` | Monotone regression backend |
| `--calc-q-from-fdr` | off | Derive q‑values from running FDRs (needs decoys) |
| `--calc-q-from-pep` | off | Derive q‑values **after** estimating PEPs |
| `--output FILE\|DIR` | *required* | Write a **target‑only** list; if DIR, a default name is used |
| `--verbose` | off | Echo all parsed parameters |

### Additional flags – **q2pep only**

| flag | default | description |
|------|---------|---------|
| `--qcol COL` | `q-value` | Column containing the input q‑values |
| `--ip` | off | Smooth PAVA step‑function with a monotone spline |
| `--ip-algo {ispline,pchip}` | `ispline` | Interpolator used with `--ip` |
| `--center-method {mean,median}` | `mean` | x‑coordinate of each PAVA block center when interpolating |

---

## Examples

### PyPI

#### q2pep

**Example 1: a target input file with q‑values**
```bash
pyisopep q2pep --target-file example/peptide.target.txt --qcol q-value --calc-q-from-pep --output example/results
```
**Example 2: separate target and decoy input files, derive q‑values from FDR first**
```bash
pyisopep q2pep --target-file example/peptide.target.txt --decoy-file example/peptide.decoy.txt --score-col score --label-col type --target-label 0 --decoy-label 1 --calc-q-from-fdr --calc-q-from-pep --output results/
```

**Example 3: a concatenated target and decoy input file with q‑values**
```bash
pyisopep q2pep --cat-file example/peptide.cat.txt --qcol q-value --score-col score --label-col type --target-label 0 --decoy-label 1 --calc-q-from-pep --output results/
```

#### d2pep

**Example 1: separate target and decoy input files**
```bash
pyisopep d2pep --target-file example/peptide.target.txt --decoy-file example/peptide.decoy.txt --score-col score --label-col type --target-label 0 --decoy-label 1 --calc-q-from-fdr --calc-q-from-pep --output results/
```

**Example 2: a concatenated target and decoy input file**
```bash
pyisopep d2pep --cat-file example/peptide.cat.txt --qcol q-value --score-col score --label-col type --target-label 0 --decoy-label 1 --calc-q-from-pep --output results/
```

### Docker image

#### Pull the Docker image from GitHub Container Registry
```bash
podman pull ghcr.io/statisticalbiotechnology/pyisotonicpep:main
```
#### q2pep
```bash
podman run --rm -it -v .:/data pyisotonicpep:main q2pep --target-file /example/peptide.target.txt --decoy-file /example/peptide.decoy.txt --score-col score --label-col type --target-label 0 --decoy-label 1 --calc-q-from-fdr --calc-q-from-pep --output /data
```

#### d2pep
```bash
podman run --rm -it -v .:/data pyisotonicpep:main d2pep --cat-file /example/peptide.cat.txt --qcol q-value --score-col score --label-col type --target-label 0 --decoy-label 1 --calc-q-from-pep --output /data
```

---

## Links
* **PyPI package:** <https://pypi.org/project/pyIsoPEP/>
* **Docker image:** <https://github.com/statisticalbiotechnology/smooth_q_to_pep/pkgs/container/pyisotonicpep>
* **Github repository:** <https://github.com/statisticalbiotechnology/smooth_q_to_pep>

---