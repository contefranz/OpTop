#' @title OpTop Package Overview
#'
#' @description
#' **OpTop** is an R package that implements the parametric approach of
#' Lewis and Grossetti (2022) for identifying the optimal number of topics of
#' a Latent Dirichlet Allocation (LDA) model, together with a family of
#' regression-style goodness-of-fit indices. Both build on the multinomial
#' structure that LDA itself assumes, so model choice and model assessment
#' rest on statistical tests rather than heuristics such as perplexity.
#'
#' The package operates on
#' **[topicmodels](https://cran.r-project.org/package=topicmodels)** fits
#' within the **[quanteda](https://quanteda.io/)** ecosystem, returns
#' **[data.table](https://rdatatable.gitlab.io/data.table/)** objects
#' throughout, and delegates its hot loops to compiled C++ so that large
#' vocabularies and long model grids remain tractable.
#'
#' @details
#' The toolkit is organized around three questions:
#'
#' - **Which K is optimal?**
#'   [optimal_topic()] evaluates a grid of fitted LDA models with the Test 1
#'   chi-square statistic of the paper (Equation 8) and selects the optimal
#'   topic count with one of three rules: the sequential adequacy scan
#'   (default), the global minimum of the standardized statistic (the rule of
#'   the published case study), or the deprecated pre-0.9.9 `"legacy"` rule.
#'   Because the chi-square reference is a yardstick rather than an exact
#'   null law, p-values can optionally be **calibrated** under the
#'   fitted-model null — `calibrate = "bootstrap"` for the exact parametric
#'   bootstrap on the collapsed envelope bins, `calibrate = "moment"` for the
#'   closed-form Haldane/Satterthwaite approximation — turning `alpha` into a
#'   genuine Type-I error rate (document lengths are supplied via
#'   `doc_lengths`).
#'
#' - **How much does a model explain?**
#'   [optop_index_se()], [optop_index_chisq()] and [optop_index_deviance()]
#'   compute regression-style \eqn{R^2} indices against the no-topics corpus
#'   baseline, at the document or word level, with micro and macro
#'   aggregations and an optional cross-document Z-test;
#'   [optop_index_table()] sweeps them across a whole grid of models.
#'
#' - **Are comparisons across K fair?**
#'   [optop_make_partition()] fixes a document-specific rare-word partition
#'   across the grid and [optop_make_baseline()] fixes the corpus baseline,
#'   so every model is evaluated on a common support.
#'
#' The vignette walks through the full workflow on the U.S. Presidential
#' Inaugural Address corpus, including the role of the envelope parameter `q`
#' and the calibrated p-values: `vignette("OpTop")`.
#'
#' The legacy stability diagnostics (`topic_stability()`,
#' `agg_topic_stability()`, `agg_document_stability()`,
#' `get_topic_models()`) are deprecated and scheduled for removal before
#' v1.0.0, together with `selection = "legacy"`.
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' Haldane, J. B. S. (1937). The exact value of the moments of the
#' distribution of chi-square. *Biometrika*, 29, 133--143.
#'
#' @name OpTop-package
#' @aliases OpTop
#' @docType package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
