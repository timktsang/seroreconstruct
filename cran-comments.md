# CRAN submission comments — seroreconstruct 1.1.3

## Resubmission

This is a resubmission addressing CRAN pre-test feedback on v1.1.2:

- **Fixed**: Replaced static `src/Makevars` with a `configure` script that
  conditionally adds `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` on Linux only.
  This resolves the `undefined symbol: dpotrf_` failure on the Debian
  (clang-21) pre-test environment while avoiding an Abort trap on macOS
  from double-linking.
- **Fixed**: Expanded acronyms in DESCRIPTION title/description (HAI, MCMC).
  Added `\value` tags to all `.Rd` files. Changed long-running examples from
  `\dontrun{}` to `\donttest{}`. All plot functions reset `par()` on exit.
- **Noted**: "Misspelled" words (Tsang, et, al, titer) are all correct —
  titer is the standard US English spelling, and Tsang et al. is an author
  citation.

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
- Ubuntu (GitHub Actions, `devel`)
- Windows Server (GitHub Actions, `release`)
