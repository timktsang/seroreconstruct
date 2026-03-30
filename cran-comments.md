# CRAN submission comments — seroreconstruct 1.1.5

## Resubmission

This is a resubmission addressing CRAN pre-test feedback on v1.1.3/v1.1.4:

- **Fixed**: Reduced total check time from >10 min to ~3 min. Tests now use
  shared pre-computed MCMC fixtures (one helper file) instead of ~49
  independent MCMC runs. Vignette iterations reduced from 5000 to 2000.
- **Fixed**: Added `inst/WORDLIST` for domain-specific terms flagged by the
  spell checker (HAI, hemagglutination, titer, Tsang, et, al). These are
  standard epidemiology/serology terminology and author citation.
- **Fixed**: Replaced GPL license badge URL in README.md that caused a
  connection timeout on the Windows pre-test server.

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
Users can override by setting the `RCPP_PARALLEL_NUM_THREADS` environment
variable before loading the package; the cap is only applied when the
variable is unset.

## Reverse dependencies

None — this is a new submission.

## Tested on

- macOS (local): R 4.2.0
- Ubuntu (GitHub Actions, `release`)
- Windows Server (GitHub Actions, `release`)
