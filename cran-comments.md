# CRAN submission comments — seroreconstruct 1.1.0

## R CMD check results

0 errors | 0 warnings | 1 note

### Note 1: New submission
```
N  checking CRAN incoming feasibility
   Maintainer: 'Tim Tsang <timkltsang@gmail.com>'
   New submission
```
This is a new submission; the package is not currently on CRAN.

## Thread policy

RcppParallel threads are capped to 2 in `.onLoad()` per CRAN policy.
Users can override via the `RCPP_PARALLEL_NUM_THREADS` environment variable.

## Reverse dependencies

None — this is a new submission.

## Tested on

- macOS (local): R 4.4.x
- Ubuntu 22.04 (GitHub Actions, `release`)
- Ubuntu 22.04 (GitHub Actions, `devel`)
- Windows Server 2022 (GitHub Actions)
