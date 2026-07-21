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

# Dense rare rows from either partition format, for the slow references and
# the structural assertions: format 2 stores the complement (non-rare word
# lists), so the mask is TRUE everywhere except at the listed words.
ref_rare_mask <- function(partition, J, W) {
  if (!identical(partition$format, 2L)) {
    return(partition$rare_mask)
  }
  mask <- matrix(TRUE, J, W)
  off <- partition$nonrare_offsets
  words <- partition$nonrare_words
  for (j in seq_len(J)) {
    if (off[j + 1] > off[j]) {
      mask[j, words[(off[j] + 1):off[j + 1]] + 1L] <- FALSE
    }
  }
  colnames(mask) <- partition$vocab
  mask
}

ref_theta_phi <- function(model) {
  if (inherits(model, "optop_theta_phi")) {
    return(list(theta = model$theta, phi = model$phi,
                K = ncol(model$theta)))
  }
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
                               reopt = "none", eps = 1e-12, min_null = 0) {
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

  rare_all <- ref_rare_mask(partition, J, ncol(N))
  for (j in seq_len(J)) {
    rare <- rare_all[j, ]
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
      dK <- sum((N_j[!rare] - E_f[!rare])^2 / E_f[!rare])
      dnull <- sum((N_j[!rare] - B_f[!rare])^2 / B_f[!rare])
      # Pearson min-bin inclusion rule: the collapsed bin enters both sides
      # only when the grid-wide flag of the partition holds (absent flag =
      # pre-0.13.0 partition = always included)
      min_ok <- if (is.null(partition$chisq_min_ok)) TRUE else
        partition$chisq_min_ok[j]
      if (isTRUE(min_ok)) {
        # format-2 convention: the min-bin mass is the unfloored sum,
        # floored ONCE (the pre-0.15.0 kernels summed floored elements)
        E_min <- max(sum(E_j[rare]), eps)
        B_min <- max(L[j] * sum(pi_row[rare]), eps)
        dK <- dK + (N_min - E_min)^2 / E_min
        dnull <- dnull + (N_min - B_min)^2 / B_min
      }
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
    r2_doc[j] <- if (dnull > 0) 1 - dK / dnull else NA_real_
  }

  # aggregation over J+ = {j : D_null(j) >= min_null}; min_null = 0 keeps
  # the legacy strict positivity rule
  valid <- if (min_null > 0) D_null >= min_null else D_null > 0
  r2_doc[!valid] <- NA_real_
  list(r2 = 1 - sum(D_K[valid]) / sum(D_null[valid]),
       r2_macro = mean(r2_doc[valid]),
       r2_doc = r2_doc,
       D_K = D_K, D_null = D_null,
       n_null_excluded = sum(!valid & is.finite(D_null)),
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
      # fitted side: Poisson unit deviance 2 * [N log(N/E) - (N - E)], with
      # the linear correction on the UNfloored expectation; for the null side
      # the correction vanishes identically (sum_j B_jw = sum_j N_jw)
      idx <- N_w > 0
      log_term <- if (any(idx)) {
        2 * sum(N_w[idx] * (log(N_w[idx]) - log(pmax(E_w, eps)[idx])))
      } else {
        0
      }
      d_w[w] <- log_term - 2 * sum(N_w - E_w)
      if (any(idx)) {
        d_w_null[w] <- 2 * sum(N_w[idx] * (log(N_w[idx]) - log(pmax(B_w, eps)[idx])))
      }
    }
  }

  # aggregation over the non-degenerate words V+ only
  valid <- d_w_null > 0
  r2_word <- rep(NA_real_, W)
  r2_word[valid] <- 1 - d_w[valid] / d_w_null[valid]
  names(r2_word) <- vocab

  list(r2_word = r2_word,
       r2_micro_word = sum((d_w_null[valid] / sum(d_w_null[valid])) *
                             r2_word[valid]),
       r2_macro_word = mean(r2_word[valid]),
       d_w = d_w, d_w_null = d_w_null,
       K = tp$K, metric = metric)
}

# Naive LEGACY reference for optimal_topic(selection = "legacy"). Mirrors the
# pre-0.9.9 semantics verbatim, one explicit loop per document:
# X_j = gamma_j %*% exp(beta), sorted descending; the envelope is truncated at
# the first position whose cumulative mass, rounded half-away-from-zero to 4
# decimals (std::round), exceeds q, with the crossing word collapsed into the
# tail; chi_j = icut * (sum((o - e)^2 / e) + pooled term). Per model the
# statistic is sum(chi_j) / sum(icut_j) and the p-value is the LOWER tail of a
# chi-square with 1 df.
#
# `W_prop` is a plain dense matrix of row-wise word proportions whose rows and
# columns are aligned with the models, exactly as the C++ core assumes.
ref_optimal_topic_legacy <- function(models, W_prop, q = 0.80) {
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
    c(topic = ncol(dtw), OpTop = stat, df = sum(icut),
      pval = stats::pchisq(stat, df = 1, lower.tail = TRUE))
  })
  as.data.frame(do.call(rbind, out))
}

# Naive reference for the calibrated Test 1 of Lewis & Grossetti (2022),
# Eq. (8): per document, the P_j important words are the smallest
# descending-sorted head whose cumulative mass strictly exceeds q (the
# crossing word is KEPT, so the collapsed tail keeps mass < 1 - q, footnote
# 5); the document contributes (P_j + 1) * Pearson over the P_j + 1 bins with
# P_j degrees of freedom (bins - 1; if the tail is empty the min bin is
# dropped, the multiplier is P_j and the df is P_j - 1). The corpus statistic
# is the RAW sum over documents, chi-square with df = sum_j P_j, and the
# p-value is the UPPER tail. The reported OpTop is the standardized statistic
# raw / df, the quantity plotted in the paper's Figure 2.
ref_optimal_topic <- function(models, W_prop, q = 0.95) {
  tps <- lapply(models, function(m) list(theta = m@gamma, phi = exp(m@beta)))
  ref_optimal_topic_tp(tps, W_prop, q = q)
}

# Same Eq. (8) reference on plain (theta, phi) pairs, for engines that do not
# store gamma/beta S4 slots (seededlda fits, NLPstudio-style wrappers). Rows
# of theta are assumed aligned with the rows of W_prop.
ref_optimal_topic_tp <- function(tps, W_prop, q = 0.95) {
  out <- lapply(tps, function(tp) {
    dtw <- tp$theta                # J x K
    tww <- t(tp$phi)               # W x K
    J <- nrow(W_prop)

    chi <- numeric(J)
    dfs <- numeric(J)
    for (j in seq_len(J)) {
      X <- drop(tww %*% dtw[j, ])
      o <- W_prop[j, ]

      ord <- order(X, decreasing = TRUE)
      X <- X[ord]
      o <- o[ord]

      cum <- cumsum(X)
      above <- which(cum > q)
      p_j <- if (length(above)) above[1] else length(X)

      head_idx <- seq_len(p_j)
      pearson <- sum((o[head_idx] - X[head_idx])^2 / X[head_idx])
      n_bins <- p_j
      if (p_j < length(X)) {
        o_tail <- sum(o[-head_idx])
        X_tail <- sum(X[-head_idx])
        pearson <- pearson + (o_tail - X_tail)^2 / X_tail
        n_bins <- p_j + 1L
      }

      chi[j] <- n_bins * pearson
      dfs[j] <- n_bins - 1L
    }

    raw <- sum(chi)
    df <- sum(dfs)
    c(topic = ncol(dtw), OpTop = raw / df, df = df,
      pval = stats::pchisq(raw, df = df, lower.tail = FALSE))
  })
  as.data.frame(do.call(rbind, out))
}

# Naive reference for the per-document envelope exported by the C++ core:
# the bin probabilities (kept-word fitted probabilities plus the collapsed
# min-bin mass) that the calibration layer consumes. Mirrors the Eq. (8)
# cutoff of ref_optimal_topic().
ref_envelope <- function(model, q = 0.95) {
  dtw <- model@gamma
  tww <- t(exp(model@beta))
  lapply(seq_len(nrow(dtw)), function(j) {
    X <- sort(drop(tww %*% dtw[j, ]), decreasing = TRUE)
    cum <- cumsum(X)
    above <- which(cum > q)
    p_j <- if (length(above)) above[1] else length(X)
    if (p_j < length(X)) {
      c(X[seq_len(p_j)], sum(X[-seq_len(p_j)]))
    } else {
      X
    }
  })
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

# Pure-R oracle for the bootstrap null: the pre-0.12.0 implementation, one
# rmultinom() call per document under the R session RNG. The compiled core
# uses its own RNG stream, so the two never match draw by draw — they must
# match in distribution, which is what the parallel tests assert.
ref_boot_null <- function(probs, doc_lengths, n_boot) {
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

# Naive reference for the held-out chain of Section 3.7 (topicmodels fits):
# explicit loops for the alignment, the OOV convention, the fold-in via
# posterior(), the held-out harmonized partition over K0 = K U {null} with
# the TRAINING baseline, the Pearson min-bin rule on the evaluation side,
# the three discrepancies, and the Macro/Micro/gap statistics with their
# standard errors.
ref_holdout <- function(models, dtm_eval, baseline, c = 1,
                        metrics = c("se", "chisq", "deviance"),
                        conf = 0.95, eps = 1e-12, min_null = 0) {
  vocab <- names(baseline$pi_glob)
  pi_tr <- unname(baseline$pi_glob)
  W <- length(vocab)
  N_raw <- as.matrix(dtm_eval)

  # alignment to the training support + out-of-support counts
  N_al <- matrix(0, nrow(N_raw), W,
                 dimnames = list(rownames(N_raw), vocab))
  hit <- colnames(N_raw)[colnames(N_raw) %in% vocab]
  N_al[, hit] <- N_raw[, hit]
  oov <- rowSums(N_raw[, setdiff(colnames(N_raw), vocab), drop = FALSE])

  keep <- rowSums(N_al) > 0
  N_al <- N_al[keep, , drop = FALSE]
  oov <- oov[keep]
  J <- nrow(N_al)
  L <- rowSums(N_al) + oov            # L_j keeps every token
  tau <- c / pmax(L, 1)

  # fold-in and fitted probabilities per model
  thetas <- lapply(models, function(m) {
    stm <- slam::as.simple_triplet_matrix(N_al)
    topicmodels::posterior(m, newdata = stm)$topics
  })
  phis <- lapply(models, function(m) ref_theta_phi(m)$phi)
  P <- lapply(seq_along(models), function(i) thetas[[i]] %*% phis[[i]])

  # harmonized union over K0 = K U {null}: min(pi_tr(w), min_K p) < tau_j;
  # the OOV pseudo-word (zero probability everywhere) is always rare
  p_min <- Reduce(pmin, P)
  p_min <- pmin(p_min, matrix(pi_tr, J, W, byrow = TRUE))
  rare <- p_min < tau

  # Pearson min-bin rule on the evaluation side (the OOV word adds no mass)
  E_min_by_model <- sapply(P, function(pm) L * rowSums(pm * rare))
  B_min <- L * as.numeric(rare %*% pi_tr)
  min_ok <- pmin(apply(E_min_by_model, 1, min), B_min) >= c

  z <- stats::qnorm(1 - (1 - conf) / 2)
  out <- list()
  for (metric in metrics) {
    per_k <- lapply(seq_along(models), function(i) {
      D_K <- numeric(J); D_0 <- numeric(J)
      for (j in seq_len(J)) {
        r_j <- rare[j, ]
        N_j <- N_al[j, ]
        E_j <- L[j] * P[[i]][j, ]
        B_j <- L[j] * pi_tr
        N_min <- sum(N_j[r_j]) + oov[j]   # OOV tokens enter the min bin

        if (metric == "se") {
          D_K[j] <- sum((N_j[!r_j] - E_j[!r_j])^2) +
            (N_min - sum(E_j[r_j]))^2
          D_0[j] <- sum((N_j[!r_j] - B_j[!r_j])^2) +
            (N_min - sum(B_j[r_j]))^2
        } else if (metric == "chisq") {
          E_f <- pmax(E_j, eps); B_f <- pmax(B_j, eps)
          dK <- sum((N_j[!r_j] - E_f[!r_j])^2 / E_f[!r_j])
          d0 <- sum((N_j[!r_j] - B_f[!r_j])^2 / B_f[!r_j])
          if (min_ok[j]) {
            # format-2 convention: unfloored sum, floored once
            E_m <- max(sum(E_j[r_j]), eps)
            B_m <- max(L[j] * sum(pi_tr[r_j]), eps)
            dK <- dK + (N_min - E_m)^2 / E_m
            d0 <- d0 + (N_min - B_m)^2 / B_m
          }
          D_K[j] <- dK; D_0[j] <- d0
        } else {
          E_m <- max(sum(E_j[r_j]), eps)
          B_m <- max(L[j] * sum(pi_tr[r_j]), eps)
          dK_min <- if (N_min == 0) 0 else 2 * N_min * log(N_min / E_m)
          d0_min <- if (N_min == 0) 0 else 2 * N_min * log(N_min / B_m)
          D_K[j] <- sum(ref_dev_contrib(N_j[!r_j], E_j[!r_j], eps)) + dK_min
          D_0[j] <- sum(ref_dev_contrib(N_j[!r_j], B_j[!r_j], eps)) + d0_min
        }
      }
      # J+ = {j : D_0(j) >= min_null}; min_null = 0 keeps strict positivity
      valid <- if (min_null > 0) D_0 >= min_null else D_0 > 0
      r <- rep(NA_real_, J)
      r[valid] <- 1 - D_K[valid] / D_0[valid]
      rv <- r[valid]; A <- D_K[valid]; B <- D_0[valid]; n <- sum(valid)
      macro <- mean(rv); macro_se <- stats::sd(rv) / sqrt(n)
      micro <- 1 - mean(A) / mean(B)
      g2 <- c(-1 / mean(B), mean(A) / mean(B)^2)
      micro_se <- sqrt(drop(t(g2) %*% stats::cov(cbind(A, B)) %*% g2) / n)
      g3 <- c(g2, -1)
      gap_se <- sqrt(drop(t(g3) %*% stats::cov(cbind(A, B, rv)) %*% g3) / n)
      list(r = r, macro = macro, macro_se = macro_se,
           ci = macro + c(-1, 1) * z * macro_se,
           micro = micro, micro_se = micro_se,
           gap = micro - macro, gap_se = gap_se, n = n,
           n_excluded = sum(!valid & is.finite(D_0)),
           D_K = D_K, D_0 = D_0)
    })
    out[[metric]] <- per_k
  }
  out
}
