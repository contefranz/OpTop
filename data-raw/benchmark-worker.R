# Worker for data-raw/benchmark-efficiency.R: one (scale, mode, threads)
# configuration per process, so GNU time can attribute a peak RSS to each.
# Not intended to be run by hand.
#
# Arguments: <fits_rds> <scale> <mode> <threads> <out_rds>
#   mode "stat"   — statistic only, compiled cores, n_threads = <threads>
#   mode "boot"   — bootstrap calibration (B = 200), n_threads = <threads>
#   mode "boot_r" — the pre-0.12.0 baseline: statistic core single-threaded
#                   plus the R rmultinom() bootstrap loop

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
} else {
  stop("unknown mode: ", mode)
}

saveRDS(list(scale = scale, mode = mode, threads = threads,
             elapsed = elapsed, result = as.data.frame(res)),
        out_rds)
