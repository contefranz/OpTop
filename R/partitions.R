#' Harmonized Rare-word Partition for OpTop indices
#'
#' Compute the fixed, document-specific rare-word sets \eqn{C_j^*} used to
#' evaluate the OpTop goodness-of-fit indices on a **common support** across a
#' grid of topic models. Following Lewis and Grossetti (2026), the union that
#' defines the harmonized set includes the no-topics baseline: a word \eqn{w}
#' is rare in document \eqn{j} if
#' \eqn{\min\bigl(\hat\pi_{\mathrm{glob}}(w),\, \min_K p^{K}_{jw}\bigr) < \tau_j,}{min(pi_glob(w), min_K p^K_jw) < tau_j,}
#' where \eqn{p^{K}_{jw}} is the fitted probability of word \eqn{w} under the
#' \eqn{K}-topic model, \eqn{\hat\pi_{\mathrm{glob}}}{pi_glob} is the corpus
#' word distribution of [optop_make_baseline()], and \eqn{\tau_j = c / L_j}
#' with \eqn{L_j} the document length.
#'
#' @param models A list of fitted topic models spanning a grid of \eqn{K}.
#'   Supports every class handled by the internal adapters:
#'   \code{topicmodels::LDA} (VEM or Gibbs) and \code{topicmodels::CTM}
#'   fits, the seededlda models and NLPstudio \code{nlp_topic_fit} objects.
#' @param dtm A document–term matrix of **counts** with rows = documents and
#'   columns = vocabulary (recommended class \code{Matrix::dgCMatrix}), or
#'   an [optop_corpus()] of count shards streamed one at a time. The rows,
#'   concatenated in shard order, must align with the models' documents.
#' @param c Positive scalar controlling the per-document threshold
#'   \eqn{\tau_j = c / L_j}. Larger values classify more words as rare.
#'   Default: \code{1}, the value the paper adopts and recommends when the
#'   deviance index is the primary measure; \code{c = 5} remains appropriate
#'   when the Pearson index is of primary interest, and sensitivity to
#'   \code{c} should be reported.
#' @param n_threads Integer; number of OpenMP threads used by the compiled
#'   kernels (default \code{1L}). Results are identical for any value; only
#'   wall time changes.
#' @param pi_glob Optional numeric vector of baseline word probabilities used
#'   as the null member of the harmonized union and in the Pearson inclusion
#'   rule. Default \code{NULL} computes the corpus distribution of
#'   \code{dtm}, the in-sample construction. Supply the \emph{training}
#'   baseline (\code{optop_make_baseline()$pi_glob}) when partitioning an
#'   evaluation corpus for held-out scoring, where the null model must come
#'   from the training sample.
#'
#' @return A list (the sparse partition, `format = 2L`) with:
#' \itemize{
#'   \item \code{nonrare_offsets}, \code{nonrare_words}: the complement of
#'   \eqn{C_j^*} in ragged form. Document \eqn{j}'s non-rare words are the
#'   \code{nonrare_words} entries (0-based word indices, ascending) between
#'   offsets \code{j} and \code{j + 1}; every other word is rare. The total
#'   size is bounded by \eqn{\sum_j L_j / c}{sum_j L_j / c}, so the
#'   structure scales with the token count rather than with \eqn{J \times W}.
#'   \item \code{vocab}: the vocabulary the word indices refer to.
#'   \item \code{L}: numeric vector of length \code{J} with document lengths
#'   \eqn{L_j = \sum_w N_{jw}}.
#'   \item \code{chisq_min_ok}: logical vector of length \code{J};
#'   \code{TRUE} when the collapsed min-bin of document \eqn{j} satisfies the
#'   Pearson inclusion rule (see Details).
#'   \item \code{chisq_min_report}: list with the number and share of
#'   documents whose min-bin is excluded and the mean observed probability
#'   mass excluded with it.
#'   \item \code{c}: the threshold constant used.
#'   \item \code{format}: the partition format tag (\code{2L}). Partitions
#'   saved by OpTop 0.13.0 to 0.14.3 (dense \code{rare_mask}) are upgraded
#'   automatically on first use, with an alert; recompute to avoid the
#'   conversion.
#' }
#'
#' @details
#' OpTop indices (SE, Pearson \eqn{\chi^2}, and Deviance) are computed on the
#' fixed support \eqn{\{w \not\in C_j^*\} \cup \{\min\}} for each document,
#' where all words in \eqn{C_j^*} are collapsed into a single "\code{min}"
#' bin. In the paper's notation the union runs over the augmented grid
#' \eqn{\mathcal{K}_0 = \mathcal{K} \cup \{\mathrm{null}\}}{K0 = K U {null}},
#' which treats the no-topics baseline as an additional model with fitted
#' probabilities \eqn{p^{\mathrm{null}}_{jw} = \hat\pi_{\mathrm{glob}}(w)}{p^null_jw = pi_glob(w)}.
#' Using a harmonized partition ensures comparability across \eqn{K} and
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
#' The construction is candidate-based and exact. Because the null baseline
#' belongs to the augmented union, a cell can be non-rare only if
#' \eqn{\hat\pi_{\mathrm{glob}}(w) \ge \tau_j}{pi_glob(w) >= tau_j}; sorting
#' the vocabulary once by \eqn{\hat\pi_{\mathrm{glob}}}{pi_glob} descending
#' makes each document's candidate set a prefix of at most \eqn{L_j / c}
#' words, and only candidates are checked against the models. Total cost is
#' \eqn{O((\sum_j L_j / c) \cdot \sum_K K)}{O((tokens / c) * sum_K K)}
#' instead of the \eqn{O(J W \sum_K K)}{O(J * W * sum_K K)} dense products
#' of earlier versions, and no \eqn{J \times W}{J x W} object of any kind is
#' materialized. The min-bin masses of the Pearson rule are complements of
#' compensated non-rare sums, exact for documents whose support collapses
#' entirely into the min-bin and accurate to
#' \eqn{O(L_j\,\epsilon_{mach})}{O(L_j * eps_machine)} otherwise. For
#' large corpora, keep \code{dtm} sparse.
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
#' str(part$nonrare_offsets)
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
optop_make_partition <- function(models, dtm, c = 1,
                                 n_threads = 1L, pi_glob = NULL) {
  stopifnot(length(models) >= 1)
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || c <= 0) {
    stop("c must be a single positive number")
  }
  if (!is.numeric(n_threads) || length(n_threads) != 1L ||
      !is.finite(n_threads) || n_threads < 1) {
    stop("n_threads must be a single integer >= 1")
  }
  n_threads <- as.integer(n_threads)
  corpus <- .optop_as_corpus_internal(dtm)
  n_sh <- corpus$n_shards
  n_models <- length(models)

  # lazy models: loader elements are materialized on demand and released
  # after each pass; in-memory fits are adapted once and cached
  is_loader <- vapply(models, is.function, logical(1))
  tp_cache <- vector("list", n_models)
  vocab <- NULL
  tp_get <- function(i) {
    if (!is.null(tp_cache[[i]])) {
      return(tp_cache[[i]])
    }
    tp <- optop_as_theta_phi(.optop_materialize_model(models[[i]]))
    if (is.null(vocab)) {
      vocab <<- colnames(tp$phi)
    } else if (!identical(colnames(tp$phi), vocab)) {
      stop(paste("all models must share the same vocabulary in the same",
                 "order"))
    }
    if (!is_loader[i]) {
      tp_cache[[i]] <<- tp
    }
    tp
  }
  tp1 <- tp_get(1L)
  W <- ncol(tp1$phi)
  J_models <- nrow(tp1$theta)

  # shard scan: document lengths and, when needed, the pooled baseline
  shard_rows <- integer(n_sh)
  L_list <- vector("list", n_sh)
  cs <- NULL
  for (s in seq_len(n_sh)) {
    m <- .optop_corpus_shard(corpus, s)
    if (ncol(m) != W) {
      stop("dtm dimensions do not match the models; align the dtm first")
    }
    shard_rows[s] <- nrow(m)
    L_list[[s]] <- Matrix::rowSums(m)
    if (is.null(pi_glob)) {
      cs_s <- Matrix::colSums(m)
      cs <- if (is.null(cs)) cs_s else cs + cs_s
    }
  }
  L <- unlist(L_list)                                 # L_j
  J <- length(L)
  if (J != J_models) {
    stop("dtm dimensions do not match the models; align the dtm first")
  }
  tau <- c / pmax(L, 1)                               # tau_j = c / L_j

  # the null baseline enters the harmonized union (K0 = K U {null}): rare iff
  # min(pi_glob(w), min_K p^K_jw) < tau_j. In sample the baseline is the
  # corpus distribution of dtm; for held-out partitions the caller supplies
  # the training baseline instead.
  if (is.null(pi_glob)) {
    pi_glob <- as.numeric(cs) / sum(cs)
  } else {
    pi_glob <- as.numeric(pi_glob)
    if (length(pi_glob) != W) {
      stop("pi_glob must have one entry per feature of dtm")
    }
  }

  # sparse construction (format 2): the partition stores the COMPLEMENT of
  # C*_j, the per-document non-rare word lists, bounded in total by
  # sum_j L_j / c (at most L_j / c words satisfy p_jw >= tau_j = c / L_j).
  # Candidates: since the null belongs to the augmented union K0, a cell can
  # be non-rare only if pi_glob(w) >= tau_j, so sorting the vocabulary by
  # pi_glob descending makes each document's candidate set a prefix
  ord <- order(pi_glob, decreasing = TRUE)
  pi_sorted <- pi_glob[ord]
  order0 <- as.integer(ord - 1L)
  prefix_len <- optop_partition_candidates_core(pi_sorted, as.numeric(tau),
                                                n_threads)
  cand_off <- c(0, cumsum(as.numeric(prefix_len)))
  keep <- as.raw(rep(1L, cand_off[J + 1]))

  # one pass per model over the document slices of each shard: a candidate
  # survives while p^K_jw >= tau_j under every model of the grid. Loaders
  # are materialized once per pass (twice per model in total).
  shard_off <- c(0, cumsum(shard_rows))
  for (i in seq_len(n_models)) {
    tp <- tp_get(i)
    if (nrow(tp$theta) != J) {
      stop("dtm dimensions do not match the models; align the dtm first")
    }
    for (s in seq_len(n_sh)) {
      th_s <- if (n_sh == 1L) tp$theta else
        tp$theta[shard_off[s] + seq_len(shard_rows[s]), , drop = FALSE]
      optop_partition_pass_core(th_s, tp$phi, order0, cand_off, keep,
                                as.numeric(tau), shard_off[s], n_threads)
    }
  }
  nr <- optop_partition_compact_core(cand_off, keep, order0, n_threads)

  # Pearson min-bin inclusion rule, decided once for the whole grid: keep
  # the collapsed bin iff min(min_K E^K_min, B_min) >= c, with the min-bin
  # masses by compensated complement, E^K_min = L * (1 - sum_NR p^K) and
  # B_min = L * (1 - sum_NR pi)
  psum_max <- rep(-Inf, J)
  for (i in seq_len(n_models)) {
    tp <- tp_get(i)
    for (s in seq_len(n_sh)) {
      idx <- shard_off[s] + seq_len(shard_rows[s])
      th_s <- if (n_sh == 1L) tp$theta else tp$theta[idx, , drop = FALSE]
      s_val <- optop_partition_sums_core(th_s, tp$phi, nr$offsets, nr$words,
                                         shard_off[s], n_threads)
      psum_max[idx] <- pmax(psum_max[idx], s_val)
    }
  }
  E_min_min <- L * pmax(0, 1 - psum_max)
  pisum <- optop_partition_pisum_core(nr$offsets, nr$words, pi_glob,
                                      n_threads)
  B_min <- L * pmax(0, 1 - pisum)
  chisq_min_ok <- pmin(E_min_min, B_min) >= c

  # a document without rare words has no min bin: the rule evaluates FALSE
  # there (so the kernel adds no spurious bin), but nothing is excluded and
  # the report counts only documents whose existing min bin is dropped
  nr_count <- diff(nr$offsets)
  has_min <- nr_count < W
  excluded <- has_min & !chisq_min_ok
  excluded_mass <- NA_real_
  if (any(excluded)) {
    obs_nr <- numeric(J)
    for (s in seq_len(n_sh)) {
      idx <- shard_off[s] + seq_len(shard_rows[s])
      Nt <- Matrix::t(.optop_corpus_shard(corpus, s))
      obs_nr[idx] <- optop_partition_obsmass_core(Nt@p, Nt@i, Nt@x,
                                                  nr$offsets, nr$words,
                                                  shard_off[s], n_threads)
    }
    N_min <- pmax(0, L - obs_nr)
    excluded_mass <- mean(N_min[excluded] / pmax(L[excluded], 1))
  }

  list(nonrare_offsets = nr$offsets,
       nonrare_words = nr$words,
       vocab = if (is.null(vocab)) .optop_corpus_vocab(corpus) else vocab,
       L = L,
       chisq_min_ok = chisq_min_ok,
       chisq_min_report = list(n_excluded = sum(excluded),
                               share = mean(excluded),
                               excluded_mass = excluded_mass),
       c = c,
       format = 2L)
}

# Slice a format-2 partition to the documents [offset + 1, offset + n]: the
# per-shard view the index engine feeds to the document kernel. The word
# list is copied for the slice (bounded by the shard's tokens / c); offsets
# are rebased to zero.
.optop_partition_slice <- function(partition, offset, n) {
  j <- offset + seq_len(n)
  off <- partition$nonrare_offsets[c(j, offset + n + 1)]
  base <- off[1]
  words <- if (off[n + 1] > base) {
    partition$nonrare_words[(base + 1):off[n + 1]]
  } else {
    integer(0)
  }
  list(nonrare_offsets = off - base,
       nonrare_words = words,
       vocab = partition$vocab,
       L = partition$L[j],
       chisq_min_ok = partition$chisq_min_ok[j],
       chisq_min_report = partition$chisq_min_report,
       c = partition$c,
       format = 2L)
}

#' Upgrade a dense-mask partition to the sparse format
#'
#' Partitions built by OpTop 0.13.0 to 0.14.3 carry the dense logical
#' `rare_mask` (J x W). The sparse format stores the complement, the
#' per-document non-rare word lists, which the format-2 kernels consume.
#' The conversion is exact; an alert recommends recomputation because the
#' dense object itself is the scale bottleneck.
#'
#' @param partition A partition from [optop_make_partition()], any format.
#'
#' @return A format-2 partition list.
#'
#' @keywords internal
.optop_partition_upgrade <- function(partition) {
  if (identical(partition$format, 2L)) {
    return(partition)
  }
  if (is.null(partition$rare_mask)) {
    stop(paste("the partition has neither the sparse fields nor a",
               "rare_mask; recompute optop_make_partition()"))
  }
  cli::cli_alert_info(paste(
    "upgrading a dense-mask partition (OpTop < 0.15.0) to the sparse",
    "format; recompute optop_make_partition() to avoid the conversion"
  ))
  mask <- partition$rare_mask
  J <- nrow(mask)
  nr_list <- vector("list", J)
  for (j in seq_len(J)) {
    nr_list[[j]] <- which(!mask[j, ]) - 1L
  }
  counts <- lengths(nr_list)
  partition$nonrare_offsets <- c(0, cumsum(as.numeric(counts)))
  partition$nonrare_words <- as.integer(unlist(nr_list, use.names = FALSE))
  partition$vocab <- colnames(mask)
  partition$rare_mask <- NULL
  partition$format <- 2L
  partition
}

# Reconstruct one document's dense rare row from the sparse partition: the
# deprecated SE re-optimization path is the only remaining consumer of a
# dense mask and stays a small-corpus code path.
.optop_partition_rare_row <- function(partition, j, W) {
  out <- rep(TRUE, W)
  o0 <- partition$nonrare_offsets[j]
  o1 <- partition$nonrare_offsets[j + 1]
  if (o1 > o0) {
    out[partition$nonrare_words[(o0 + 1):o1] + 1L] <- FALSE
  }
  out
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
#'   columns = vocabulary (recommended class \code{Matrix::dgCMatrix}), or
#'   an [optop_corpus()] of count shards, whose counts are pooled.
#' @param smooth_lambda Nonnegative additive smoothing constant:
#'   \eqn{\hat\pi^{\lambda}_{\mathrm{glob}}(w) \propto N_{\cdot w} + \lambda}{pi^lambda_glob(w) proportional to N_.w + lambda}.
#'   Default \code{0} (no smoothing, the in-sample construction). A positive
#'   value keeps every baseline probability strictly positive, the smoothing
#'   option of the paper's held-out support convention; sensitivity to
#'   \eqn{\lambda}{lambda} should be reported.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{pi_glob}: numeric vector of length \code{W} with
#'   \eqn{\pi_{\mathrm{glob}}(w) = N_{\cdot w} / N_{\cdot \cdot}}, i.e., the
#'   corpus-level term probabilities (smoothed when
#'   \code{smooth_lambda > 0}).
#' }
#'
#' @details
#' The returned \code{pi_glob} is used by the exported index functions to form
#' per-document baseline counts on the same evaluation support defined by
#' \code{optop_make_partition()}. Keeping this baseline fixed across \eqn{K}
#' enables interpretable, regression-style \eqn{R^2} measures: the null and
#' the saturated models are the two interpretative boundaries between which
#' every fitted \eqn{K}-topic model lies. For held-out evaluation, compute
#' the baseline on the \emph{training} corpus and pass it to
#' [optop_index_holdout()]; words absent from the training sample have zero
#' baseline probability unless \code{smooth_lambda > 0}.
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
optop_make_baseline <- function(dtm, smooth_lambda = 0) {
  if (!is.numeric(smooth_lambda) || length(smooth_lambda) != 1L ||
      !is.finite(smooth_lambda) || smooth_lambda < 0) {
    stop("smooth_lambda must be a single nonnegative number")
  }
  # pi^lambda_glob(w) = (N_.w + lambda) / sum_w (N_.w + lambda), with the
  # counts pooled across the shards when dtm is an optop_corpus
  if (.optop_is_corpus(dtm)) {
    cs <- NULL
    for (s in seq_len(dtm$n_shards)) {
      cs_s <- Matrix::colSums(.optop_corpus_shard(dtm, s))
      cs <- if (is.null(cs)) cs_s else cs + cs_s
    }
    N_tot <- cs + smooth_lambda
    vocab <- .optop_corpus_vocab(dtm)
  } else {
    N_tot <- Matrix::colSums(dtm) + smooth_lambda
    vocab <- colnames(dtm)
  }
  pi_glob <- as.numeric(N_tot) / sum(N_tot)
  names(pi_glob) <- vocab   # <-- add names for safe alignment
  list(pi_glob = pi_glob)
}
