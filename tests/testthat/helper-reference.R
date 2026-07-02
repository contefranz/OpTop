# Naive reference implementation of the OpTop discrepancy indices.
#
# This mirrors the definitions of the methodological paper as directly as
# possible: explicit per-document (or per-word) loops on the harmonized support
# {w not in C*_j} U {min}, with no blocking, no rowSums decomposition tricks and
# no vectorized shortcuts. The optimized package code must reproduce these
# numbers exactly (up to floating-point noise); any disagreement is a red flag
# to be reported, not papered over.
#
# Conventions replicated from the package code:
# * chisq: expected counts are floored at eps element-wise BEFORE any use;
#   the collapsed "min" bin uses the sum of the floored elements, itself
#   floored at eps. The baseline min bin, however, is L_j * sum(pi[rare])
#   computed from the UNfloored pi, then floored at eps.
# * deviance: contributions are 2 * N * log(N / max(E, eps)) with the
#   0 * log(0) = 0 convention; the min bin uses the sum of UNfloored expected
#   counts, floored at eps, and log(max(N_min, eps) / E_min).
# * se: no eps flooring anywhere.

ref_theta_phi <- function(model) {
  p <- topicmodels::posterior(model)
  list(theta = p$topics, phi = p$terms, K = ncol(p$topics))
}

# Expected counts E (J x W): E_jw = L_j * sum_k theta_jk phi_kw
ref_expected <- function(theta, phi, L) {
  E <- theta %*% phi
  E * L  # column-major recycling scales row j by L[j]
}

ref_dev_contrib <- function(N, E, eps = 1e-12) {
  out <- 2 * N * log(N / pmax(E, eps))
  out[N == 0] <- 0
  out
}

# Document-level indices, one slow loop per document.
ref_index_document <- function(model, dtm, partition, baseline,
                               metric = c("se", "chisq", "deviance"),
                               reopt = "none", eps = 1e-12) {
  metric <- match.arg(metric)
  tp <- ref_theta_phi(model)
  vocab <- colnames(tp$phi)
  pi_row <- baseline$pi_glob[vocab]
  N <- as.matrix(dtm)
  L <- partition$L
  J <- nrow(N)
  E <- ref_expected(tp$theta, tp$phi, L)

  D_K <- numeric(J)
  D_null <- numeric(J)
  r2_doc <- numeric(J)

  for (j in seq_len(J)) {
    rare <- partition$rare_mask[j, ]
    N_j <- N[j, ]
    E_j <- E[j, ]
    B_j <- L[j] * pi_row
    N_min <- sum(N_j[rare])

    if (metric == "se") {
      if (identical(reopt, "se")) {
        # lambda blending of E towards B on the harmonized support
        num <- sum((N_j[!rare] - B_j[!rare]) * (E_j[!rare] - B_j[!rare])) +
          (N_min - sum(B_j[rare])) * (sum(E_j[rare]) - sum(B_j[rare]))
        den <- sum((E_j[!rare] - B_j[!rare])^2) +
          (sum(E_j[rare]) - sum(B_j[rare]))^2
        lambda <- if (den <= 0) 0 else max(0, min(1, num / den))
        E_j <- lambda * E_j + (1 - lambda) * B_j
      }
      dK <- sum((N_j[!rare] - E_j[!rare])^2) + (N_min - sum(E_j[rare]))^2
      dnull <- sum((N_j[!rare] - B_j[!rare])^2) + (N_min - sum(B_j[rare]))^2
    } else if (metric == "chisq") {
      E_f <- pmax(E_j, eps)
      B_f <- pmax(B_j, eps)
      E_min <- max(sum(E_f[rare]), eps)
      B_min <- max(L[j] * sum(pi_row[rare]), eps)
      dK <- sum((N_j[!rare] - E_f[!rare])^2 / E_f[!rare]) +
        (N_min - E_min)^2 / E_min
      dnull <- sum((N_j[!rare] - B_f[!rare])^2 / B_f[!rare]) +
        (N_min - B_min)^2 / B_min
    } else { # deviance
      E_min <- max(sum(E_j[rare]), eps)
      B_min <- max(L[j] * sum(pi_row[rare]), eps)
      dK_min <- if (N_min == 0) 0 else 2 * N_min * log(max(N_min, eps) / E_min)
      dnull_min <- if (N_min == 0) 0 else 2 * N_min * log(max(N_min, eps) / B_min)
      dK <- sum(ref_dev_contrib(N_j[!rare], E_j[!rare], eps)) + dK_min
      dnull <- sum(ref_dev_contrib(N_j[!rare], B_j[!rare], eps)) + dnull_min
    }

    D_K[j] <- dK
    D_null[j] <- dnull
    r2_doc[j] <- if (dnull > 0) 1 - dK / dnull else 0
  }

  list(r2 = 1 - sum(D_K) / sum(D_null),
       r2_macro = mean(r2_doc[D_null > 0]),
       r2_doc = r2_doc,
       D_K = D_K, D_null = D_null,
       K = tp$K, metric = metric)
}

# Word-level indices, one slow loop per word. The harmonized partition plays
# no role here beyond providing document lengths L.
ref_index_word <- function(model, dtm, partition, baseline,
                           metric = c("se", "chisq", "deviance"),
                           eps = 1e-12) {
  metric <- match.arg(metric)
  tp <- ref_theta_phi(model)
  vocab <- colnames(tp$phi)
  pi_row <- baseline$pi_glob[vocab]
  N <- as.matrix(dtm)
  L <- partition$L
  W <- ncol(N)
  E <- ref_expected(tp$theta, tp$phi, L)

  d_w <- numeric(W)
  d_w_null <- numeric(W)

  for (w in seq_len(W)) {
    N_w <- N[, w]
    E_w <- E[, w]
    B_w <- L * pi_row[w]
    if (metric == "se") {
      d_w[w] <- sum((N_w - E_w)^2)
      d_w_null[w] <- sum((N_w - B_w)^2)
    } else if (metric == "chisq") {
      E_f <- pmax(E_w, eps)
      B_f <- pmax(B_w, eps)
      d_w[w] <- sum((N_w - E_f)^2 / E_f)
      d_w_null[w] <- sum((N_w - B_f)^2 / B_f)
    } else { # deviance
      idx <- N_w > 0
      if (any(idx)) {
        d_w[w] <- 2 * sum(N_w[idx] * (log(N_w[idx]) - log(pmax(E_w, eps)[idx])))
        d_w_null[w] <- 2 * sum(N_w[idx] * (log(N_w[idx]) - log(pmax(B_w, eps)[idx])))
      }
    }
  }

  r2_word <- ifelse(d_w_null > 0, 1 - d_w / d_w_null, 0)
  names(r2_word) <- vocab

  list(r2_word = r2_word,
       r2_micro_word = sum((d_w_null / sum(d_w_null)) * r2_word),
       r2_macro_word = mean(r2_word[d_w_null > 0]),
       K = tp$K, metric = metric)
}

# Naive reference for optimal_topic(). Mirrors the semantics of
# src/optimal_topic_core.cpp verbatim, one explicit loop per document:
# X_j = gamma_j %*% exp(beta), sorted descending; the envelope is truncated at
# the first position whose cumulative mass, rounded half-away-from-zero to 4
# decimals (std::round), exceeds q; the truncated tail is pooled into a single
# bin; chi_j = icut * (sum((o - e)^2 / e) + pooled term). Per model the
# statistic is sum(chi_j) / sum(icut_j) and the p-value is the LOWER tail of a
# chi-square with 1 df (a documented oddity preserved on purpose — see
# AUDIT.md).
#
# `W_prop` is a plain dense matrix of row-wise word proportions whose rows and
# columns are aligned with the models, exactly as the C++ core assumes.
ref_optimal_topic <- function(models, W_prop, q = 0.80) {
  out <- lapply(models, function(m) {
    dtw <- m@gamma                 # J x K
    tww <- t(exp(m@beta))          # W x K
    J <- nrow(W_prop)

    chi <- numeric(J)
    icut <- numeric(J)
    for (j in seq_len(J)) {
      X <- drop(tww %*% dtw[j, ])
      o <- W_prop[j, ]

      ord <- order(X, decreasing = TRUE)
      X <- X[ord]
      o <- o[ord]

      cum <- cumsum(X)
      # floor(v * 1e4 + 0.5) reproduces std::round() for the positive values
      # cum takes; R's round() would use banker's rounding at exact halves.
      above <- which(floor(cum * 1e4 + 0.5) / 1e4 > q)
      ic <- if (length(above)) above[1] - 1L else length(X)

      head_idx <- seq_len(ic)
      o_tail <- sum(o[-head_idx])
      X_tail <- sum(X[-head_idx])
      chi[j] <- ic * (sum((o[head_idx] - X[head_idx])^2 / X[head_idx]) +
                        (o_tail - X_tail)^2 / X_tail)
      icut[j] <- ic
    }

    stat <- sum(chi) / sum(icut)
    c(topic = ncol(dtw), OpTop = stat,
      pval = stats::pchisq(stat, df = 1, lower.tail = TRUE))
  })
  as.data.frame(do.call(rbind, out))
}

# Reference for the cross-document Z-test.
ref_ztest <- function(r2_doc_valid) {
  J <- length(r2_doc_valid)
  m <- mean(r2_doc_valid)
  se <- sqrt(sum((r2_doc_valid - m)^2) / (J - 1)) / sqrt(J)
  z <- m / se
  list(z = z,
       pval = stats::pnorm(z, lower.tail = FALSE),
       se = se,
       ci = m + c(-1, 1) * stats::qnorm(0.975) * se,
       J = J)
}
