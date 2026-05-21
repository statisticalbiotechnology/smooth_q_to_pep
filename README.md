# `pyIsoPEP`: Isotonic PEP Estimator

## Overview
`pyIsoPEP` provides both a [**Docker image**](https://github.com/statisticalbiotechnology/smooth_q_to_pep/pkgs/container/pyisotonicpep) (Docker-, Podman-, and Apptainer-compatible) and a [**Python API**](https://pypi.org/project/pyIsoPEP/) with a unified interface for estimating smooth, non‑decreasing Posterior Error Probabilities (PEPs) directly from empirical null models using isotonic regression for target identifications in shotgun proteomics.

Six estimation combinations are available — three input methods × two regressors — mirroring the same design used in [**Percolator**](http://percolator.ms/) (C++):

| # | method | x-axis | input | regressor | CLI |
|---|--------|--------|-------|-----------|-----|
| 1 | **q2pep** | rank | target q-values | I-Spline | `q2pep` |
| 2 | **q2pep** | rank | target q-values | PAVA | `q2pep --pava` |
| 3 | **qns2pep** | score | target q-values + scores | I-Spline | `q2pep --score-based` |
| 4 | **qns2pep** | score | target q-values + scores | PAVA | `q2pep --score-based --pava` |
| 5 | **tdc2pep** | score | TDC targets + decoys | I-Spline | `d2pep` |
| 6 | **tdc2pep** | score | TDC targets + decoys | PAVA | `d2pep --pava` |

**q2pep / qns2pep** fit PEPs from target q-values only; only target statistics are reported.  
**tdc2pep** fits `p(decoy | score)` over all PSMs and converts to PEP via `p / (1 − p)`; only target PEPs are reported.  
**qns2pep** and **tdc2pep** place I-Spline knots at quantiles of the *score* distribution, concentrating spline flexibility where observations are densest.

Optional post-processing derives q-values from estimated PEPs for both DDA and DIA data.

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
# or
apptainer build pyisopep.sif docker://ghcr.io/statisticalbiotechnology/pyisotonicpep:main
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
```

Workflow‑specific help:

```bash
pyisopep q2pep -h
pyisopep d2pep -h
```

### Common flags (valid for *both* subcommands)

| flag | default | description |
|------|---------|---------|
| `--cat-file FILE`, `--target-file FILE` | *required* | Select concatenated vs. separate input |
| `--decoy-file FILE` | — | Decoy file (required in separate input mode) |
| `--score-col COL` | `score` | Score column |
| `--label-col COL` | `label` | Column distinguishing targets/decoys in a concatenated file |
| `--target-label STR` | `target` | String marking target rows |
| `--decoy-label STR` | `decoy` | String marking decoy rows |
| `--pava` | off | Use PAVA regression instead of the default I-Spline |
| `--calc-q-from-fdr` | off | Derive q‑values from running FDRs (needs decoys) |
| `--calc-q-from-pep` | off | Derive q‑values **after** estimating PEPs |
| `-k`, `--keep-input-order` | off | Write output in original input order |
| `--output FILE\|DIR` | *required* | Write a **target‑only** list; if DIR, a default name is used |
| `--verbose` | off | Echo all parsed parameters |

### Additional flags – **q2pep only**

| flag | default | description |
|------|---------|---------|
| `--qcol COL` | `q-value` | Column containing the input q‑values |
| `--score-based` | off | Use raw scores as independent variable instead of rank (selects **qns2pep**); requires `--score-col` and decoy information |
| `--trim-plateaus` | off | Enable trimming of leading/trailing q‑value plateaus before isotonic regression (plateau trimming is off by default) |
| `--no-pseudo-count` | off | Disable the Jeffreys-style smoothing (value `0.5`, distributed as `0.5 / n_mid` per rank) that is added to raw PEPs before isotonic regression by default. Smoothing improves log-ratio calibration at very small q. |

---

## Examples

### PyPI

#### Combination 1 — q2pep, I-Spline (rank-based, default)

```bash
pyisopep q2pep --target-file example/peptide.target.txt --qcol q-value --calc-q-from-pep --output example/results
```

#### Combination 2 — q2pep, PAVA (rank-based)

```bash
pyisopep q2pep --target-file example/peptide.target.txt --qcol q-value --pava --calc-q-from-pep --output results/
```

#### Combination 3 — qns2pep, I-Spline (score-based q-values)

```bash
pyisopep q2pep --score-based \
  --target-file example/peptide.target.txt \
  --decoy-file  example/peptide.decoy.txt \
  --qcol q-value --score-col score \
  --calc-q-from-pep --output results/
```

#### Combination 4 — qns2pep, PAVA

```bash
pyisopep q2pep --score-based --pava \
  --target-file example/peptide.target.txt \
  --decoy-file  example/peptide.decoy.txt \
  --qcol q-value --score-col score \
  --calc-q-from-pep --output results/
```

#### Combination 5 — tdc2pep, I-Spline (score-based TDC)

```bash
pyisopep d2pep \
  --target-file example/peptide.target.txt \
  --decoy-file  example/peptide.decoy.txt \
  --score-col score --calc-q-from-pep --output results/
```

#### Combination 6 — tdc2pep, PAVA

```bash
pyisopep d2pep --pava \
  --target-file example/peptide.target.txt \
  --decoy-file  example/peptide.decoy.txt \
  --score-col score --calc-q-from-pep --output results/
```

### Docker image

```bash
podman pull ghcr.io/statisticalbiotechnology/pyisotonicpep:main
```

```bash
# q2pep (rank-based, default)
podman run --rm -it -v .:/data pyisotonicpep:main \
  q2pep --target-file /example/peptide.target.txt \
  --decoy-file /example/peptide.decoy.txt \
  --score-col score --qcol q-value \
  --calc-q-from-fdr --calc-q-from-pep --output /data

# tdc2pep (score-based TDC)
podman run --rm -it -v .:/data pyisotonicpep:main \
  d2pep \
  --cat-file /example/peptide.cat.txt \
  --score-col score --label-col type --target-label 0 --decoy-label 1 \
  --calc-q-from-pep --output /data
```

---

## Python API

```python
import numpy as np
from pyIsoPEP.IsotonicPEP import IsotonicPEP

iso = IsotonicPEP()          # I-Spline by default
iso_pava = IsotonicPEP(pava=True)   # PAVA

# --- Combinations 1 & 2: q2pep (rank-based) ---
_, _, pep, q2 = iso.pep_regression(
    q_values=q_sorted, method="q2pep",
    pava=False,   # True for PAVA
    calc_q_from_pep=True,
)

# --- Combinations 3 & 4: qns2pep (score-based q-values) ---
_, _, pep, q2 = iso.pep_regression(
    q_values=q_sorted, obs=obs_all,   # obs shape (n, 2): [score, label]
    method="qns2pep",
    pava=False,   # True for PAVA
    calc_q_from_pep=True,
)

# --- Combinations 5 & 6: tdc2pep (score-based TDC) ---
_, _, pep, q2 = iso.pep_regression(
    obs=obs_all,
    method="tdc2pep",
    pava=False,   # True for PAVA
    calc_q_from_pep=True,
)

# Direct access — all PSMs (targets + decoys):
df_obs = iso.process_obs(obs_all)
pep_all = iso.tdc_to_pep(df_obs)            # I-Spline
pep_all = iso.tdc_to_pep(df_obs, pava=True) # PAVA

# Unified q_to_pep — rank-based (scores=None) or score-based (scores provided):
pep_rank  = iso.q_to_pep(q_sorted)
pep_score = iso.q_to_pep(q_sorted, scores=score_sorted)
```

---

## Links
* **PyPI package:** <https://pypi.org/project/pyIsoPEP/>
* **Docker image:** <https://github.com/statisticalbiotechnology/smooth_q_to_pep/pkgs/container/pyisotonicpep>
* **Github repository:** <https://github.com/statisticalbiotechnology/smooth_q_to_pep>

---
