[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/contefranz/OpTop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/contefranz/OpTop/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/contefranz/OpTop/graph/badge.svg)](https://app.codecov.io/gh/contefranz/OpTop)
[![release](https://img.shields.io/badge/release-v0.15.0-blue.svg)](https://github.com/contefranz/OpTop/releases/tag/v0.15.0)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)

# OpTop

OpTop provides principled, statistically grounded tools to select the optimal number of topics and 
to assess goodness-of-fit for probabilistic topic models. 
It implements fast parametric tests and discrepancy indices that make comparisons across 
different topic counts directly comparable. The methodology is developed on the multinomial 
structure of the Latent Dirichlet Allocation (LDA) family and applies to any fitted model 
that returns document-level word probabilities: supported engines include **topicmodels** 
(LDA with VEM or Gibbs, and CTM), **seededlda**, and **NLPstudio**, which can also be mixed 
within one model grid.

### What It Does

- **Optimal K selection**  
  Fast, parametric tests based on a Pearson-type statistic identify the topic count that best describes the corpus (`optimal_topic()`), with optional bootstrap or moment-matched calibration of the p-values under the fitted-model null (v0.10.0). The pre-0.9.8 redundant-topic diagnostics were removed in v0.16.0 after a long deprecation.

- **Model goodness-of-fit**  
  Regression-style indices for topic models (e.g., SE, Pearson-$\chi^2$, and deviance), summarized 
  as micro/macro $R^2$ analogues (`optop_index_se()`, `optop_index_chisq()`, `optop_index_deviance()`), 
  plus a convenience table across a grid (`optop_index_table()`).

- **Held-out validation and residual diagnostics**  
  As of v0.14.0, `optop_index_holdout()` scores independent evaluation documents with 
  confidence intervals for the average held-out fit, `optop_gain_table()` selects the 
  smallest sufficient K from paired adjacent gains, and `optop_moment_test()` tests 
  held-out residuals for structured bias along frequency or fit-based vocabulary 
  strata.

- **Harmonized comparisons across K**  
  Rare words are collapsed via a fixed, document-specific partition so all indices are evaluated on a common support (`optop_make_partition()`), with a fixed corpus baseline for fair normalization (`optop_make_baseline()`).

- **Performance at scale**  
  As of v0.15.0 the package evaluates corpora of tens of millions of documents with
  bounded memory: the corpus crosses the compiled boundary as zero-copy sparse views,
  the harmonized partition stores per-document non-rare word lists whose size scales
  with the token count rather than with documents times features, corpora beyond the
  2^31 - 1 nonzero container cap of a single sparse matrix are supplied as document
  shards through `optop_corpus()` (streamed one at a time, with per-document results
  bit-identical to an unsharded run), and model grids too large for memory are
  supplied as loader functions materialized one model at a time.
  Core routines are implemented in C/C++ for large vocabularies and model grids. As of
  v0.9.7, the `optimal_topic()` core works on blocked BLAS matrix products and pairs
  documents with the fitted models by identifier, so the dfm row order no longer matters.
  As of v0.12.0, both the Test 1 statistic and the bootstrap calibration are parallelized
  with OpenMP, controlled by the `n_threads` argument of `optimal_topic()`. Results are
  bit-identical for any thread count and reproducible under a seed; on toolchains without
  OpenMP support the routines run single-threaded. The bootstrap selects its sampling
  algorithm per document (alias-table token draws on wide envelopes, conditional-binomial
  draws otherwise) and envelope ties are resolved deterministically, so results are
  identical across platforms, including under the tied fitted probabilities that Gibbs
  estimators produce. On a 4-core machine the bootstrap-calibrated pass runs up to 45
  times faster than in v0.11.0 (a corpus of 10,000 documents and 20,000 features drops
  from 1 h 47 m to 142 s); see the vignette section "A Note On Computational Efficiency"
  for the full simulation study.

- **Current support and extensions**  
  As of v0.11.0, `optimal_topic()` and the discrepancy indices accept, through the 
  internal adapter family (`optop_as_theta_phi()`), **topicmodels** fits (LDA with 
  VEM or Gibbs, and `CTM`), all three **seededlda** models (`textmodel_lda()`, 
  `textmodel_seededlda()`, `textmodel_seqlda()`), and **NLPstudio** fits 
  (`nlp_topic_fit`). Engines can be mixed within one grid fitted on the same corpus: 
  documents and features are aligned by identifier, per model. Supporting a new 
  engine requires a single adapter method that returns the fitted document-topic 
  and topic-word probabilities; *WarpLDA* fits from **text2vec** are supported 
  through the `optop_warplda()` wrapper.
  
### Authors

- [Francesco Grossetti](https://contefranz.github.io/) 

  _Assistant Professor of Accounting Analytics and Data Science_  
  Department of Accounting | Bocconi University  
  Fellow at Bocconi Institute for Data Science and Analytics ([BIDSA](https://www.bidsa.unibocconi.eu/wps/wcm/connect/Site/Bidsa/Home))  
  Contact: francesco.grossetti@unibocconi.it

- [Craig M. Lewis](https://business.vanderbilt.edu/bio/craig-lewis/)

  _Madison S. Wigginton Professor of Finance, Emeritus_  
  Owen Graduate School of Management | Vanderbilt University  
  Contact: craig.lewis@vanderbilt.edu  
  
### References

1. Lewis, C. M. and Grossetti, F. (2022). 
[*A statistical approach for optimal topic model identification*](https://www.jmlr.org/papers/volume23/19-297/19-297.pdf). 
Journal of Machine Learning Research, 23(58), 1–20.
2. Grün, B. and Hornik, K. (2011). 
[*topicmodels: An R package for fitting topic models.*](https://www.jstatsoft.org/article/view/v040i13) 
Journal of Statistical Software, 40, 1-30.
3. Blei, D. M., Ng, A. Y., and Jordan, M. I. (2003). 
[*Latent Dirichlet Allocation*.](https://www.jmlr.org/papers/volume3/blei03a/blei03a.pdf)
Journal of Machine Learning Research, 3(Jan):993–1022.
4. Grossetti F. (2026). [*NLPstudio: Tools for Social Science Text Analysis, Corpus Management, and Topic Modeling*](https://contefranz.github.io/NLPstudio). R package version
  1.1.1.
  

### Installation

OpTop compiles C++ code (via Rcpp/RcppArmadillo) and links your R installation's
BLAS/LAPACK, so installing from GitHub builds from source and needs a C/C++ toolchain
**and gfortran** in place first.

- **Windows.** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching
  your R version (Rtools43 for R 4.3, Rtools44 for R 4.4, and so on). It provides the
  C/C++ compiler, gfortran, `make`, and OpenMP; nothing else is required and
  multithreading works out of the box.
- **macOS.** Install the Command Line Tools for Xcode with `xcode-select --install` (C/C++
  compiler and `make`), then the official R gfortran toolchain from
  [mac.r-project.org/tools](https://mac.r-project.org/tools/) (for example
  `gfortran-14.2-universal.pkg`). gfortran is required to link the BLAS/LAPACK Fortran
  runtime; without it the build fails at the link step with a missing `-lgfortran`. Apple's
  default clang has no OpenMP, so OpTop runs single-threaded on a stock toolchain (results
  are identical, only serial); to enable multithreading install LLVM/`libomp` and point
  `~/.R/Makevars` at it.
- **Linux.** Install a C/C++ compiler, gfortran, and the R headers, for example on
  Debian/Ubuntu `sudo apt-get install build-essential gfortran r-base-dev`. GCC provides
  OpenMP automatically.

```r
# From GitHub
# install.packages("pak")
pak::pkg_install("contefranz/OpTop")
```

The package vignette presents the complete workflow (optimal-K selection, the role of `q`,
the discrepancy indices, and a simulation with known truth): `vignette("OpTop")`.

### Quick Start

```r
library(OpTop)
library(quanteda)
library(topicmodels)

# 0) Get a quanteda corpus object
corpus_texts = data_corpus_inaugural

# 1) Build a weighted dfm (row-wise proportions) for optimal-K selection
dfm_counts   <- dfm(tokens(corpus_texts))                    # counts
weighted_dfm <- dfm_weight(dfm_counts, scheme = "prop")      # proportions

# 2) Fit a grid of topic models (e.g., K = 10:50 by 5); LDA via VEM is shown,
#    and any adapter-supported engine works (topicmodels, seededlda, NLPstudio)
set.seed(123)
K_grid    <- seq(10, 50, by = 5)
VEM_models <- lapply(
  K_grid, function(k) {
    lda = topicmodels::LDA(x = dfm_counts, k = k, method = "VEM")
  }
) 

# 3) Choose the optimal K (selection: "sequential" adequacy scan by default,
#    "min" for the paper's global-minimum rule; verbose = TRUE reports progress)
res_opt   <- optimal_topic(topic_models = VEM_models, 
                           weighted_dfm = weighted_dfm,
                           q = 0.95, 
                           alpha = 0.05, 
                           selection = "sequential",
                           do_plot = TRUE,
                           verbose = TRUE)

# 4) Goodness-of-fit across K (use counts here)
part      <- optop_make_partition(models = VEM_models, dtm = dfm_counts, c = 1)
base      <- optop_make_baseline(dtm = dfm_counts)
tab       <- optop_index_table(models = VEM_models, 
                               dtm = dfm_counts,
                               metrics = c("se","chisq","deviance"),
                               partition = part, 
                               baseline = base, 
                               macro = TRUE)
tab
```

### At very large scale (optional)

Nothing above changes for corpora that fit in a single dfm. When a corpus outgrows one
sparse matrix (the container caps at 2^31 - 1 nonzero entries) or the model grid outgrows
memory, two optional inputs take over; everything else, and every result, stays the same:

```r
# corpora beyond a single dfm: document shards, streamed one at a time
corp <- optop_corpus(shard_files, reader = function(p) qs2::qs_read(p))

# grids beyond memory: loader functions, materialized one model at a time
loaders <- lapply(fit_files, function(f) { force(f); function() qs2::qs_read(f) })

res  <- optimal_topic(loaders, corp_weighted, n_threads = 20)
part <- optop_make_partition(loaders, corp, c = 1)
```

Sharded evaluation is exact, not an approximation: per-document results are bit-identical
to an unsharded run and a one-shard corpus reproduces the plain call bit for bit,
calibrated p-values included. See `?optop_corpus` for the full contract.

### Reproducibility across platforms

OpTop itself is deterministic: envelope ties are broken in a fixed order, the bootstrap
draws from its own internal random-number stream, and all reductions run in a fixed
order, so for a given set of fitted models it returns identical results on Windows,
macOS, and Linux and for any `n_threads`. The small differences you may see when running
the example on different machines come from the model-fitting step upstream:
`topicmodels::LDA(method = "VEM")` is a numerical optimization whose result depends on the
platform's BLAS/LAPACK and math library, so `set.seed()` fixes the initialization but not
the last-digit rounding, and the fitted probabilities differ slightly across platforms.
To compare results exactly across machines, fit once and share the fitted model objects
(OpTop returns identical output for identical inputs), or use the same BLAS on both.

### Bug Reporting

Bugs and issues can be reported at
[https://github.com/contefranz/OpTop/issues](https://github.com/contefranz/OpTop/issues).


***
  
