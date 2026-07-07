#' Harmonized Rare-word Partition for OpTop indices
#'
#' Compute the fixed, document-specific rare-word sets \eqn{C_j^*} used to
#' evaluate the OpTop goodness-of-fit indices on a **common support** across a
#' grid of topic models. Following Lewis and Grossetti (2026), the union that
#' defines the harmonized set includes the no-topics baseline: a word \eqn{w}
#' is rare in document \eqn{j} if
#' \eqn{\min\bigl(\hat\pi_{\mathrm{glob}}(w),\, \min_K i^{K}_{jw}\bigr) < \tau_j,}{min(pi_glob(w), min_K i^K_jw) < tau_j,}
#' where \eqn{i^{K}_{jw}} is the fitted probability of word \eqn{w} under the
#' \eqn{K}-topic model, \eqn{\hat\pi_{\mathrm{glob}}}{pi_glob} is the corpus
#' word distribution of [optop_make_baseline()], and \eqn{\tau_j = c / L_j}
#' with \eqn{L_j} the document length.
#'
#' @param models A list of fitted topic models spanning a grid of \eqn{K}.
#'   Supports every class handled by the internal adapters:
#'   \code{topicmodels::LDA} (VEM or Gibbs) and \code{topicmodels::CTM}
#'   fits, the seededlda models and NLPstudio \code{nlp_topic_fit} objects.
#' @param dtm A document–term matrix of **counts** with rows = documents and
#'   columns = vocabulary. Recommended class is \code{Matrix::dgCMatrix}.
#' @param c Positive scalar controlling the per-document threshold
#'   \eqn{\tau_j = c / L_j}. Larger values classify more words as rare.
#'   Default: \code{1}, the value the paper adopts and recommends when the
#'   deviance index is the primary measure; \code{c = 5} remains appropriate
#'   when the Pearson index is of primary interest, and sensitivity to
#'   \code{c} should be reported.
#' @param block Integer; number of terms to process per block when multiplying
#'   \eqn{\Theta \Phi} to control memory usage for large vocabularies.
#'   Default: \code{5000}.
#' @param n_threads Integer; number of OpenMP threads used by the compiled
#'   kernels (default \code{1L}). Results are identical for any value; only
#'   wall time changes.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{rare_mask}: logical matrix \code{J x W}; \code{TRUE} where
#'   word \code{w} is rare in document \code{j} (belongs to \eqn{C_j^*}).
#'   \item \code{L}: numeric vector of length \code{J} with document lengths
#'   \eqn{L_j = \sum_w N_{jw}}.
#'   \item \code{chisq_min_ok}: logical vector of length \code{J};
#'   \code{TRUE} when the collapsed min-bin of document \eqn{j} satisfies the
#'   Pearson inclusion rule (see Details).
#'   \item \code{chisq_min_report}: list with the number and share of
#'   documents whose min-bin is excluded and the mean observed probability
#'   mass excluded with it.
#'   \item \code{c}: the threshold constant used.
#' }
#'
#' @details
#' OpTop indices (SE, Pearson \eqn{\chi^2}, and Deviance) are computed on the
#' fixed support \eqn{\{w \not\in C_j^*\} \cup \{\min\}} for each document,
#' where all words in \eqn{C_j^*} are collapsed into a single "\code{min}"
#' bin. Using a harmonized partition ensures comparability across \eqn{K} and
#' prevents support-driven artifacts in goodness-of-fit curves. Because the
#' baseline enters the union, every non-min expected count satisfies
#' \eqn{E^{K}_{jw} \ge c} and \eqn{B_{jw} \ge c} on the active support, which
#' bounds all Pearson denominators away from zero.
#'
#' **Pearson min-bin inclusion rule.** For Pearson-type summaries the
#' collapsed min-bin of document \eqn{j} is included only when
#' \eqn{\min\bigl(\min_K E^{K}_{j,\min},\, B_{j,\min}\bigr) \ge c}{min(min_K E^K_j,min, B_j,min) >= c}
#' and is excluded otherwise, from both the fitted and the null discrepancy
#' and for all models simultaneously. The decision is taken here, once for
#' the whole grid, and exposed as \code{chisq_min_ok}; the share of documents
#' affected and the excluded observed mass are stored in
#' \code{chisq_min_report} and reported by the index functions. The rule does
#' not affect the deviance and squared-error families.
#'
#' **Grid dependence.** \eqn{C_j^*} is a union over the model grid, so the
#' indices computed on this support depend on the pair (grid, \code{c}):
#' extending the grid weakly enlarges \eqn{C_j^*} and can change every
#' reported index. Fix the grid before evaluation and compare indices across
#' studies only under a common grid and \code{c}.
#'
#' @section Computational notes:
#' The compiled kernel forms
#' \eqn{\min\bigl(\hat\pi_{\mathrm{glob}}, \min_K \Theta^{(K)} \Phi^{(K)}\bigr)}{min(pi_glob, min_K Theta Phi)}
#' one vocabulary block at a time: the baseline seeds the running minimum and
#' a single BLAS product per model per block feeds the in-place update, so
#' neither the \code{J x W x |K|} tensor nor the full \code{J x W}
#' running-minimum matrix is ever materialized. A second blocked pass of the
#' same cost class accumulates the rare-set fitted mass per model for the
#' Pearson inclusion rule. For large corpora, keep \code{dtm} sparse.
#'
#' @seealso
#' \code{\link{optop_make_baseline}},
#' \code{\link{optop_index_se}},
#' \code{\link{optop_index_chisq}},
#' \code{\link{optop_index_deviance}},
#' \code{\link{optop_index_table}}
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' library(topicmodels)
#'
#' # toy sparse DTM (J=50, W=1000)
#' set.seed(1)
#' ij <- cbind(
#'   sample(50, 5000, TRUE),
#'   sample(1000, 5000, TRUE)
#' )
#' x  <- sample(1:3, 5000, TRUE)
#' dtm <- sparseMatrix(i = ij[,1], j = ij[,2], x = x, dims = c(50, 1000))
#'
#' # fit a small grid of LDA models
#' m5  <- LDA(dtm, k = 5,  method = "VEM", control = list(seed = 42))
#' m10 <- LDA(dtm, k = 10, method = "VEM", control = list(seed = 42))
#'
#' part <- optop_make_partition(list(m5, m10), dtm, c = 1)
#' str(part$rare_mask)
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' @export
optop_make_partition <- function(models, dtm, c = 1, block = 5000,
                                 n_threads = 1L) {
  stopifnot(length(models) >= 1)
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || c <= 0) {
    stop("c must be a single positive number")
  }
  if (!is.numeric(n_threads) || length(n_threads) != 1L ||
      !is.finite(n_threads) || n_threads < 1) {
    stop("n_threads must be a single integer >= 1")
  }
  L <- Matrix::rowSums(dtm)                           # L_j
  tau <- c / pmax(L, 1)                               # τ_j
  theta_phi <- lapply(models, optop_as_theta_phi)
  J <- nrow(theta_phi[[1]]$theta); W <- ncol(theta_phi[[1]]$phi)
  if (nrow(dtm) != J || ncol(dtm) != W) {
    stop("dtm dimensions do not match the models; align the dtm first")
  }

  # the corpus baseline enters the harmonized union: rare iff
  # min(pi_glob(w), min_K i^K_jw) < tau_j
  N_tot <- Matrix::colSums(dtm)
  pi_glob <- as.numeric(N_tot) / sum(N_tot)

  thetas <- lapply(theta_phi, `[[`, "theta")
  phis <- lapply(theta_phi, `[[`, "phi")
  rare_mask <- matrix(NA, nrow = J, ncol = W,
                      dimnames = list(rownames(theta_phi[[1]]$theta),
                                      colnames(theta_phi[[1]]$phi)))
  optop_partition_fill_core(rare_mask, thetas, phis,
                            pi_glob,
                            as.numeric(tau),
                            as.integer(block),
                            as.integer(n_threads))

  # Pearson min-bin inclusion rule, decided once for the whole grid: keep
  # the collapsed bin iff min(min_K E^K_min, B_min) >= c
  rare_i_sum <- optop_partition_minmass_core(rare_mask, thetas, phis,
                                             as.integer(block),
                                             as.integer(n_threads))
  E_min_min <- L * apply(rare_i_sum, 1, min)
  B_min <- L * as.numeric(rare_mask %*% pi_glob)
  chisq_min_ok <- pmin(E_min_min, B_min) >= c

  # a document without rare words has no min bin: the rule evaluates FALSE
  # there (so the kernel adds no spurious bin), but nothing is excluded and
  # the report counts only documents whose existing min bin is dropped
  has_min <- rowSums(rare_mask) > 0
  excluded <- has_min & !chisq_min_ok
  excluded_mass <- NA_real_
  if (any(excluded)) {
    trip <- Matrix::summary(methods::as(dtm, "CsparseMatrix"))
    hit <- rare_mask[cbind(trip$i, trip$j)]
    N_min <- numeric(J)
    if (any(hit)) {
      agg <- rowsum(trip$x[hit], trip$i[hit])
      N_min[as.integer(rownames(agg))] <- agg[, 1]
    }
    excluded_mass <- mean(N_min[excluded] / pmax(L[excluded], 1))
  }

  list(rare_mask = rare_mask,
       L = L,
       chisq_min_ok = chisq_min_ok,
       chisq_min_report = list(n_excluded = sum(excluded),
                               share = mean(excluded),
                               excluded_mass = excluded_mass),
       c = c)
}

#' Global Corpus Baseline for OpTop Indices
#'
#' Compute the global word distribution \eqn{\pi_{\mathrm{glob}}} from a corpus
#' of counts. This baseline defines the "no-topics" model used to normalize
#' OpTop indices via \eqn{B_{jw} = L_j \, \pi_{\mathrm{glob}}(w)} and
#' \eqn{B_{j,\min} = L_j \sum_{w \in C_j^*} \pi_{\mathrm{glob}}(w)}. It is the
#' analogue of the intercept-only model in regression: every document follows
#' the same corpus-level word distribution, scaled by its length.
#'
#' @param dtm A document–term matrix of **counts** with rows = documents and
#'   columns = vocabulary. Recommended class is \code{Matrix::dgCMatrix}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{pi_glob}: numeric vector of length \code{W} with
#'   \eqn{\pi_{\mathrm{glob}}(w) = N_{\cdot w} / N_{\cdot \cdot}}, i.e., the
#'   corpus-level term probabilities.
#' }
#'
#' @details
#' The returned \code{pi_glob} is used by the exported index functions to form
#' per-document baseline counts on the same evaluation support defined by
#' \code{optop_make_partition()}. Keeping this baseline fixed across \eqn{K}
#' enables interpretable, regression-style \eqn{R^2} measures: the null and
#' the saturated models are the two interpretative boundaries between which
#' every fitted \eqn{K}-topic model lies.
#'
#' @seealso
#' \code{\link{optop_make_partition}},
#' \code{\link{optop_index_se}},
#' \code{\link{optop_index_chisq}},
#' \code{\link{optop_index_deviance}},
#' \code{\link{optop_index_table}}
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' set.seed(1)
#' dtm <- rsparsematrix(100, 5000, density = 0.002, rand.x = function(n) rpois(n, 2))
#' base <- optop_make_baseline(dtm)
#' sum(base$pi_glob)  # ~ 1
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' @export
optop_make_baseline <- function(dtm) {
  N_tot <- Matrix::colSums(dtm)
  L_tot <- sum(N_tot)
  pi_glob <- as.numeric(N_tot) / L_tot
  names(pi_glob) <- colnames(dtm)   # <-- add names for safe alignment
  list(pi_glob = pi_glob)
}
