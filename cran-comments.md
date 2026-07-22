# cran-comments

## Submission note

This is a new submission (first CRAN release of OpTop).

## Test environments

- local: macOS (Apple Silicon), R 4.6.1
- GitHub Actions: ubuntu-latest (R release, oldrel-1, devel),
  macos-latest (R release), windows-latest (R release), all with
  `--as-cran --run-donttest`
- GitHub Actions: ASAN/UBSAN sanitizer job on the compiled kernels
- win-builder: devel and release
- macbuilder: release

## R CMD check results

0 errors | 0 warnings | 1 note

- New submission (first release). The spell-check flag on "Grossetti" in
  the DESCRIPTION is an author surname.

## Notes for the reviewers

- The vignette rebuilds from a small precomputed results bundle shipped in
  `inst/extdata` (17 KB), so no topic models are fitted at check time and
  the vignette knits in seconds.
- Compiled OpenMP kernels default to a single thread; examples, tests, and
  the vignette never use more than 2.

## Reverse dependencies

There are no reverse dependencies; this is the first release.
