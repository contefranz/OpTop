# Null calibration for the Test 1 statistic.
#
# Both calibrations run under the conditional null "the fitted K-topic model
# generated the data" with theta and phi held fixed at their estimates: each
# document is Multinomial(N_j, I_j) over the vocabulary, with I_j the model's
# fitted word probabilities. Because the per-document envelope (the sorted
# fitted probabilities, the cutoff at q, the collapsed min bin) depends only
# on the model, and a multinomial collapsed over bins is multinomial on the
# collapsed probabilities, the null distribution of the statistic is fully
# determined by the per-document bin probabilities exported by the C++ core —
# no full-vocabulary simulation is ever needed.
#
# Performance note: the loops below touch only C-level primitives (rmultinom,
# colSums, vector arithmetic) at ~n_boot x df flops per model. Should corpora
# with J >> 10^4 documents ever make the per-document R loop the bottleneck,
# this file is the designated C++/OpenMP port — the envelope export already
# provides all of its inputs.

#' Split the flattened envelope into per-document bins
#'
#' The C++ core (`optimal_topic_core(..., return_envelope = TRUE)`) exports
#' the envelope as one concatenated numeric vector of bin probabilities plus
#' the number of bins per document; this helper splits it back into a list
#' with one bin-probability vector per document, the shape the calibration
#' routines consume.
#'
#' @param bin_probs Numeric vector; the concatenated per-document bin
#'   probabilities (kept-word fitted probabilities followed by the collapsed
#'   min-bin mass).
#' @param bin_counts Integer vector; the number of bins of each document.
#'
#' @return A list of numeric vectors, one per document, each summing to 1.
#'
#' @keywords internal
.optop_split_envelope <- function(bin_probs, bin_counts) {
  split(bin_probs, rep(seq_along(bin_counts), bin_counts))
}

#' Parametric bootstrap draw of the null Test 1 statistic
#'
#' Draws `n_boot` replicates of the statistic under the conditional
#' fitted-model null by simulating counts directly on the collapsed envelope
#' bins, \eqn{c_j \sim \mathrm{Multinomial}(N_j, p_j)}{c_j ~ Multinomial(N_j,
#' p_j)} — exact, because a multinomial collapsed over bins is multinomial on
#' the collapsed probabilities. Each replicate accumulates
#' \eqn{T^\ast = \sum_j k_j \sum_b (c_{jb}/N_j - p_{jb})^2 / p_{jb}}{T* =
#' sum_j k_j sum_b (c_jb/N_j - p_jb)^2 / p_jb}, where the multiplier
#' \eqn{k_j} is the number of bins of document \eqn{j}, matching the
#' statistic's convention (\eqn{P_j + 1} with a min bin, \eqn{P_j} without).
#' One `rmultinom()` call per document yields all replicates at once.
#'
#' @param probs List of per-document bin-probability vectors, as produced by
#'   `.optop_split_envelope()`.
#' @param doc_lengths Numeric vector of document lengths \eqn{N_j}, aligned
#'   with `probs`.
#' @param n_boot Integer; number of bootstrap replicates.
#'
#' @return A numeric vector of length `n_boot` with the null replicates
#'   \eqn{T^\ast}{T*}.
#'
#' @keywords internal
.optop_boot_null <- function(probs, doc_lengths, n_boot) {
  T_null <- numeric(n_boot)
  for (j in seq_along(probs)) {
    p <- probs[[j]]
    n <- doc_lengths[j]
    cnt <- stats::rmultinom(n_boot, size = n, prob = p)
    dev <- (cnt / n - p)^2 / p
    T_null <- T_null + length(p) * colSums(dev)
  }
  T_null
}

#' Exact null moments of the Test 1 statistic and their Satterthwaite match
#'
#' Computes the null mean and variance of the statistic in closed form and
#' the scaled chi-square \eqn{a\,\chi^2_\nu}{a * chisq(nu)} matching them.
#' Per document, with \eqn{k} bins of probabilities \eqn{p_b} and length
#' \eqn{n}, the count-based Pearson statistic \eqn{X^2} has exact multinomial
#' moments (Haldane, 1937)
#' \eqn{E[X^2] = k - 1}{E[X^2] = k - 1} and
#' \eqn{\mathrm{Var}[X^2] = 2(k - 1) + (\sum_b 1/p_b - k^2 - 2k + 2)/n}{Var[X^2] = 2(k - 1) + (sum_b 1/p_b - k^2 - 2k + 2)/n},
#' and the document's contribution to the statistic is
#' \eqn{S_j = (k_j/N_j)\,X^2_j}{S_j = (k_j/N_j) X^2_j}, so its moments scale
#' by \eqn{k_j/N_j} and \eqn{(k_j/N_j)^2}. Documents are independent, hence
#' the corpus moments are sums, and the Satterthwaite match is
#' \eqn{a = \sigma^2/(2\mu)}{a = sigma^2/(2 mu)},
#' \eqn{\nu = 2\mu^2/\sigma^2}{nu = 2 mu^2/sigma^2}; the calibrated p-value
#' is the upper tail of \eqn{T/a} on \eqn{\nu}{nu} degrees of freedom.
#'
#' @param probs List of per-document bin-probability vectors, as produced by
#'   `.optop_split_envelope()`.
#' @param doc_lengths Numeric vector of document lengths \eqn{N_j}, aligned
#'   with `probs`.
#'
#' @return A list with the null mean `mu`, the null variance `sigma2`, and
#'   the Satterthwaite scale `a` and degrees of freedom `nu`.
#'
#' @references
#' Haldane, J. B. S. (1937). The exact value of the moments of the
#' distribution of chi-square. *Biometrika*, 29, 133--143.
#'
#' @keywords internal
.optop_moment_null <- function(probs, doc_lengths) {
  mu <- 0
  sigma2 <- 0
  for (j in seq_along(probs)) {
    p <- probs[[j]]
    n <- doc_lengths[j]
    k <- length(p)
    scale_j <- k / n
    mu <- mu + scale_j * (k - 1)
    sigma2 <- sigma2 +
      scale_j^2 * (2 * (k - 1) + (sum(1 / p) - k^2 - 2 * k + 2) / n)
  }
  list(mu = mu, sigma2 = sigma2,
       a = sigma2 / (2 * mu), nu = 2 * mu^2 / sigma2)
}
