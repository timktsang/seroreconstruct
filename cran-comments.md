# CRAN submission comments — seroreconstruct 1.1.1

## Resubmission

This is a resubmission addressing CRAN pre-test feedback on v1.1.0:

- **Fixed**: Variable Length Arrays (VLA) in C++ replaced with `std::vector<double>` (3 occurrences in `mcmc_function.cpp`). This resolves the `-Wvla` warning on Windows/GCC.
- **Fixed**: Removed obsolete `CXX_STD = CXX11` from `Makevars` and `Makevars.win`, and removed `C++11` from `SystemRequirements`. The package compiles fine under the default C++ standard.
- **Noted**: "Misspelled" words (HAI, Tsang, et, al, titer) are all correct — HAI is a standard immunology abbreviation (hemagglutination inhibition), titer is the standard US English spelling, and Tsang et al. is an author citation.

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
