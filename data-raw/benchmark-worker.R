# Worker for data-raw/benchmark-efficiency.R: one (scale, mode, threads)
# configuration per process, so GNU time can attribute a peak RSS to each.
# Not intended to be run by hand.
#
# Arguments: <fits_rds> <scale> <mode> <threads> <out_rds>
#   mode "stat"      — statistic only, compiled cores, n_threads = <threads>
#   mode "boot"      — bootstrap calibration (B = 200), n_threads = <threads>
#   mode "boot_r"    — the pre-0.12.0 baseline: statistic core single-threaded
#                      plus the R rmultinom() bootstrap loop
#   mode "indices"   — optop_make_partition() + three-metric document-level
#                      optop_index_table() on reconstructed counts, compiled
#                      engine, n_threads = <threads>
#   mode "indices_r" — the pre-engine baseline: the vectorized R partition
#                      (full J x W running-minimum matrix) and the blocked R
#                      index implementation, reimplemented verbatim

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
fits_rds <- args[[1L]]
scale <- args[[2L]]
mode <- args[[3L]]
threads <- as.integer(args[[4L]])
out_rds <- args[[5L]]

suppressMessages(library(OpTop))

dat <- readRDS(fits_rds)  # the per-scale cache written by the driver
models <- dat$models
wdfm <- dat$wdfm
dl <- dat$doc_lengths
n_boot <- 200L

if (mode == "stat") {
  elapsed <- system.time(
    res <- optimal_topic(models, wdfm, n_threads = threads,
                         do_plot = FALSE, verbose = FALSE)
  )[["elapsed"]]
} else if (mode == "boot") {
  elapsed <- system.time(
    res <- optimal_topic(models, wdfm, calibrate = "bootstrap",
                         n_boot = n_boot, doc_lengths = dl, seed = 13,
                         n_threads = threads, do_plot = FALSE,
                         verbose = FALSE)
  )[["elapsed"]]
} else if (mode == "boot_r") {
  # the 0.11.0 bootstrap: envelope export per model, then one rmultinom()
  # call per document with the Pearson reduction in vectorized R
  elapsed <- system.time({
    docs <- as.character(quanteda::docid(wdfm))
    dfm_t <- Matrix::t(methods::as(wdfm, "dgCMatrix"))
    set.seed(13)
    rows <- lapply(models, function(m) {
      tp <- OpTop:::optop_as_theta_phi(m)
      doc_map <- match(docs, tp$docs) - 1L
      out <- OpTop:::optimal_topic_core(tp$theta, tp$phi, dfm_t, 0.95,
                                        doc_map, TRUE, 1L)
      probs <- OpTop:::.optop_split_envelope(out$bin_probs, out$bin_counts)
      T_obs <- out$stat[1L, 2L]
      T_null <- numeric(n_boot)
      for (j in seq_along(probs)) {
        p <- probs[[j]]
        n <- dl[j]
        cnt <- stats::rmultinom(n_boot, size = n, prob = p)
        dev <- (cnt / n - p)^2 / p
        T_null <- T_null + length(p) * colSums(dev)
      }
      data.frame(topic = out$stat[1L, 1L], OpTop = T_obs / out$stat[1L, 3L],
                 df = out$stat[1L, 3L],
                 pval = (1 + sum(T_null >= T_obs)) / (n_boot + 1))
    })
    res <- do.call(rbind, rows)
  })[["elapsed"]]
} else if (mode %in% c("indices", "indices_r")) {
  # Counts reconstructed from the cached proportion dfm: the fits cache holds
  # wdfm = counts / L, so round(prop * L) recovers the integer counts exactly.
  dtm <- methods::as(wdfm, "dgCMatrix")
  dtm@x <- round(dtm@x * dl[dtm@i + 1L])
  dtm <- Matrix::drop0(dtm)

  # The default c = 5 collapses the harmonized support entirely on these
  # flat synthetic corpora (every word of every document falls below
  # tau_j = c / L_j from medium-1 upward, see the tau degeneracy note in
  # dev/AUDIT.md), which would make the index values cancellation noise.
  # c = 0.2 keeps a non-trivial support at every scale (no fully collapsed
  # document, median rare share 0.84-0.97) so the correctness gate compares
  # values that carry signal; the arithmetic volume, and therefore the
  # timing, does not depend on c.
  c_rare <- 0.2

  if (mode == "indices") {
    t_part <- system.time(
      partition <- optop_make_partition(models, dtm, c = c_rare,
                                        n_threads = threads)
    )[["elapsed"]]
    baseline <- optop_make_baseline(dtm)
    t_table <- system.time(
      res <- optop_index_table(models, dtm,
                               metrics = c("se", "chisq", "deviance"),
                               macro = TRUE, partition = partition,
                               baseline = baseline, n_threads = threads)
    )[["elapsed"]]
    elapsed <- t_part + t_table
  } else {
    # The pre-engine implementation, reimplemented verbatim: partition via a
    # full J x W running-minimum matrix and blocked pmin(); indices via the
    # blocked dense decomposition D = D_full - D_rare + D_min, one metric at
    # a time per model, with the null side hoisted across the grid as the
    # old optop_index_table() did.
    eps <- 1e-12
    theta_phi <- lapply(models, OpTop:::optop_as_theta_phi)
    J <- nrow(theta_phi[[1L]]$theta)
    W <- ncol(theta_phi[[1L]]$phi)

    t_part <- system.time({
      L <- Matrix::rowSums(dtm)
      tau <- c_rare / pmax(L, 1)
      i_min <- matrix(Inf, nrow = J, ncol = W)
      colnames(i_min) <- colnames(theta_phi[[1L]]$phi)
      rownames(i_min) <- rownames(theta_phi[[1L]]$theta)
      for (tp in theta_phi) {
        for (start in seq(1L, W, by = 5000L)) {
          end <- min(start + 5000L - 1L, W)
          I <- tp$theta %*% tp$phi[, start:end, drop = FALSE]
          i_min[, start:end] <- pmin(i_min[, start:end], I)
        }
      }
      rare_mask <- i_min < tau
    })[["elapsed"]]
    rm(i_min)

    t_table <- system.time({
      pi_row <- as.numeric(Matrix::colSums(dtm))
      pi_row <- pi_row / sum(pi_row)
      block_docs <- max(100L, min(J, floor(5e8 / (W * 8))))
      B_min_all <- L * as.numeric(rare_mask %*% pi_row)

      # one document-block sweep: model side (theta/phi) or, when theta is
      # NULL, the baseline side B = outer(L, pi_row)
      doc_disc <- function(metric, theta, phi) {
        D <- numeric(J)
        for (start in seq(1L, J, by = block_docs)) {
          end <- min(start + block_docs - 1L, J)
          j_idx <- start:end
          N_block <- as.matrix(dtm[j_idx, , drop = FALSE])
          rare_block <- rare_mask[j_idx, , drop = FALSE]
          if (is.null(theta)) {
            E_block <- outer(L[j_idx], pi_row)
            E_min <- B_min_all[j_idx]
          } else {
            E_block <- (theta[j_idx, , drop = FALSE] %*% phi) * L[j_idx]
            E_min <- rowSums(rare_block * E_block)
          }
          N_min <- rowSums(rare_block * N_block)
          if (metric == "se") {
            contrib <- (N_block - E_block)^2
            D_min <- (N_min - E_min)^2
          } else if (metric == "chisq") {
            # the old implementation floors E element-wise first, so the
            # model-side min bin sums the floored elements; the baseline min
            # bin stays the unfloored L_j * sum(pi[rare]), floored once
            E_block <- pmax(E_block, eps)
            if (!is.null(theta)) E_min <- rowSums(rare_block * E_block)
            contrib <- (N_block - E_block)^2 / E_block
            E_min <- pmax(E_min, eps)
            D_min <- (N_min - E_min)^2 / E_min
          } else {
            contrib <- 2 * N_block * log(N_block / pmax(E_block, eps))
            contrib[N_block == 0] <- 0
            D_min <- 2 * N_min * log(pmax(N_min, eps) / pmax(E_min, eps))
            D_min[N_min == 0] <- 0
          }
          D[j_idx] <- rowSums(contrib) - rowSums(rare_block * contrib) + D_min
        }
        D
      }

      metrics <- c("se", "chisq", "deviance")
      D_null <- lapply(metrics, function(mtr) doc_disc(mtr, NULL, NULL))
      names(D_null) <- metrics
      rows <- lapply(theta_phi, function(tp) {
        # add_baseline_topic = TRUE, the old SE/chisq default: a zero-weight
        # baseline row that leaves results unchanged but rides in the BLAS
        phi_aug <- rbind(tp$phi, pi_row)
        theta_aug <- cbind(tp$theta, rep(0, J))
        res <- list(K = tp$K)
        for (mtr in metrics) {
          if (mtr == "deviance") {
            D_K <- doc_disc(mtr, tp$theta, tp$phi)
          } else {
            D_K <- doc_disc(mtr, theta_aug, phi_aug)
          }
          D0 <- D_null[[mtr]]
          r2_doc <- ifelse(D0 > 0, 1 - D_K / D0, 0)
          col <- c(se = "R2_SE", chisq = "R2_chisq", deviance = "R2_dev")[[mtr]]
          res[[col]] <- 1 - sum(D_K) / sum(D0)
          res[[paste0(col, "_macro")]] <- mean(r2_doc[D0 > 0])
        }
        as.data.frame(res, check.names = FALSE)
      })
      res <- do.call(rbind, rows)
    })[["elapsed"]]
    elapsed <- t_part + t_table
  }
} else {
  stop("unknown mode: ", mode)
}

out <- list(scale = scale, mode = mode, threads = threads,
            elapsed = elapsed, result = as.data.frame(res))
if (exists("t_part")) {
  out$elapsed_partition <- t_part
  out$elapsed_table <- t_table
}
saveRDS(out, out_rds)
