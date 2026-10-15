# ATLAS-Ridge

Reference implementation of the adaptive generalized M-ridge estimator of

> Junaid, M.; Dar, I.S.; Iftikhar, H.; Qureshi, M.; Canas Rodrigues, P.
> An Adaptive Generalized M-Ridge Estimator with Bounded Influence for
> Contaminated Multicollinear Data. *Mathematics*, 2026.

## Requirements

R (>= 4.3) with `robustbase` and `MASS`:

    install.packages(c("robustbase", "MASS"))

## Contents

| File | Purpose |
|---|---|
| `atlas_timings.R` | the estimator and the ten competitors, implemented from Sections 3 and 4 of the paper, together with the timing harness |
| `sessionInfo.txt` | machine and package versions behind the reported timings |

## Reproducing the timing table

    Rscript atlas_timings.R

About 20 minutes on the machine recorded in `sessionInfo.txt` (Intel Core
i5-6600 @ 3.30 GHz, single-threaded). The script calls `set.seed(20261015)`
before any random draw, so repeated runs give identical results.

Four blocks are printed: elapsed time per fit at three problem sizes, the
decomposition of the proposed estimator into its pipeline stages, a check
of the estimation error against the simulation study, and the session
information.

Timing uses a batched clock: each measurement times a batch of repeated
calls and divides by the batch size, so that no entry falls below the
resolution of the system clock.

## Reproducibility note

The script was run on two machines under different operating systems and
different versions of R. Estimates agreed to four decimal places and the
stage decomposition agreed to within two percentage points.
