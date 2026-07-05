[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/contefranz/OpTop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/contefranz/OpTop/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/contefranz/OpTop/graph/badge.svg)](https://app.codecov.io/gh/contefranz/OpTop)
[![release](https://img.shields.io/badge/release-v0.11.0-blue.svg)](https://github.com/contefranz/OpTop/releases/tag/v0.11.0)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)

# OpTop

OpTop provides principled, statistically grounded tools to select the optimal number of topics and 
to assess goodness-of-fit for Latent Dirichlet Allocation (LDA) models. 
It implements fast parametric tests and discrepancy indices that make comparisons across 
different topic counts directly comparable.

### What It Does

- **Optimal K selection**  
  Fast, parametric tests based on a Pearson-type statistic identify the topic count that best describes the corpus (`optimal_topic()`), with optional bootstrap or moment-matched **calibration of the p-values** under the fitted-model null (v0.10.0). The legacy redundant-topic diagnostics (`topic_stability()`, `agg_topic_stability()`, `agg_document_stability()`) are deprecated as of v0.9.8 and scheduled for removal.

- **Model goodness-of-fit**  
  Regression-style indices for topic models (e.g., SE, Pearson-$\chi^2$, and deviance), summarized 
  as micro/macro $R^2$ analogues (`optop_index_se()`, `optop_index_chisq()`, `optop_index_deviance()`), 
  plus a convenience table across a grid (`optop_index_table()`).

- **Harmonized comparisons across K**  
  Rare words are collapsed via a fixed, document-specific partition so all indices are evaluated on a common support (`optop_make_partition()`), with a fixed corpus baseline for fair normalization (`optop_make_baseline()`).

- **Performance at scale**  
  Core routines are implemented in C/C++ for large vocabularies and model grids. As of
  v0.9.7, the `optimal_topic()` core works on blocked BLAS matrix products and pairs
  documents with the fitted models by identifier, so the dfm row order no longer matters.
  As of v0.12.0, both the Test 1 statistic and the bootstrap calibration are parallelized
  with OpenMP, controlled by the `n_threads` argument of `optimal_topic()`. Results are
  bit-identical for any thread count and reproducible under a seed; on toolchains without
  OpenMP support the routines run single-threaded. On a 4-core machine the
  bootstrap-calibrated pass runs up to 5.9 times faster than in v0.11.0; see the vignette
  section "A Note On Computational Efficiency" for the full simulation study.

- **Current support and extensions**  
  As of v0.11.0, `optimal_topic()` and the discrepancy indices accept — through the 
  internal adapter family (`optop_as_theta_phi()`) — **topicmodels** fits (LDA with 
  VEM or Gibbs, and `CTM`), all three **seededlda** models (`textmodel_lda()`, 
  `textmodel_seededlda()`, `textmodel_seqlda()`), and **NLPstudio** fits 
  (`nlp_topic_fit`). Engines can be mixed within one grid fitted on the same corpus: 
  documents and features are aligned by identifier, per model. An adapter for 
  *WarpLDA* from **text2vec** is planned.
  
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
  

### Installation

```r
# From GitHub
# install.packages("remotes")
remotes::install_github("contefranz/OpTop")
```

A full walkthrough — optimal-K selection, the role of `q`, the discrepancy indices, and a
simulation with known truth — is in the package vignette: `vignette("OpTop")`.

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

# 2) Fit a grid of LDA models (e.g., K = 10:50 by 5)
set.seed(123)
K_grid    <- seq(10, 50, by = 5)
VEM_models <- lapply(
  K_grid, function(k) {
    lda = topicmodels::LDA(x = dfm_counts, k = k, method = "VEM")
  }
) 

# 3) Choose the optimal K (selection: "sequential" adequacy scan by default,
#    "min" for the paper's global-minimum rule; verbose = TRUE reports progress)
res_opt   <- optimal_topic(lda_models = VEM_models, 
                           weighted_dfm = weighted_dfm,
                           q = 0.95, 
                           alpha = 0.05, 
                           selection = "sequential",
                           do_plot = TRUE,
                           verbose = TRUE)

# 4) Goodness-of-fit across K (use counts here)
part      <- optop_make_partition(models = VEM_models, dtm = dfm_counts, c = 5)
base      <- optop_make_baseline(dtm = dfm_counts)
tab       <- optop_index_table(models = VEM_models, 
                               dtm = dfm_counts,
                               metrics = c("se","chisq","deviance"),
                               partition = part, 
                               baseline = base, 
                               macro = TRUE)
tab
```

### Bug Reporting

Bugs and issues can be reported at
[https://github.com/contefranz/OpTop/issues](https://github.com/contefranz/OpTop/issues).


***
  
