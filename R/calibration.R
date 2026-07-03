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

# Split the flattened envelope returned by the C++ core into a list of
# per-document bin-probability vectors.
.optop_split_envelope <- function(bin_probs, bin_counts) {
  split(bin_probs, rep(seq_along(bin_counts), bin_counts))
}

# Parametric bootstrap draw of the null statistic: n_boot replicates of
# T = sum_j m_j * sum_b (c_b / N_j - p_b)^2 / p_b with c ~ Multinomial(N_j, p)
# drawn directly on the collapsed bins (exact, by multinomial aggregation).
# The multiplier m_j equals the number of bins, matching the statistic's
# convention (P_j + 1 with a min bin, P_j without).
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

# Exact null moments of the statistic plus the Satterthwaite scaled
# chi-square a * chi^2_nu matching them. Per document, with k bins of
# probabilities p and length n, the count-based Pearson X^2 has exact
# multinomial moments (Haldane, 1937)
#   E[X^2]   = k - 1
#   Var[X^2] = 2 (k - 1) + (sum_b 1/p_b - k^2 - 2 k + 2) / n
# and the document's contribution is S_j = (m_j / n) * X^2 with m_j = k, so
# its moments scale by (m_j / n) and (m_j / n)^2. Documents are independent,
# hence the corpus moments are sums; a = sigma^2 / (2 mu), nu = 2 mu^2 /
# sigma^2, and the calibrated p-value is the upper tail of T / a on nu df.
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
