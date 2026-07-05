# Efficiency study for the v0.12.0 C++/OpenMP release: speed, thread
# scaling and peak memory of optimal_topic() on synthetic corpora, with the
# pre-0.12.0 R bootstrap as baseline. Feeds the "Computational efficiency"
# section of the vignette.
#
# Run from the package root, with the current OpTop installed and GNU time
# available at /usr/bin/time:
#   Rscript data-raw/benchmark-efficiency.R
#
# Design:
# * three corpus scales, VEM grids of 10 models each (topicmodels::LDA with
#   capped EM iterations: the benchmark exercises the evaluation pipeline,
#   not the fitting, and flatter fits give larger envelopes — a conservative
#   setting for the bootstrap timings);
# * each (scale, mode, threads) configuration runs in a fresh R process via
#   data-raw/benchmark-worker.R under /usr/bin/time -v, so wall time comes
#   from system.time() inside the worker and peak RSS from the OS;
# * modes: "stat" (statistic only), "boot" (bootstrap calibration, B = 200),
#   "boot_r" (the 0.11.0 R bootstrap reimplemented verbatim as baseline);
# * correctness is asserted alongside: thread sweeps must be bit-identical,
#   and the C++ bootstrap p-values must sit within Monte-Carlo resolution of
#   the R baseline.
#
# Model fits are cached per scale in data-raw/benchmark-fits-<scale>.rds
# (gitignored); delete a cache file to refit that scale. Results go to
# inst/extdata/optop-benchmark-results.rds.

suppressMessages({
  library(OpTop)
  library(quanteda)
})

scales <- list(
  small = list(J = 200L, W = 2000L),
  medium = list(J = 500L, W = 5000L),
  large = list(J = 1000L, W = 10000L)
)
K_grid <- seq(2L, 20L, by = 2L)
em_iter <- 30L
threads_sweep <- c(1L, 2L, 4L)
n_boot <- 200L

fits_rds_for <- function(scale) {
  file.path("data-raw", sprintf("benchmark-fits-%s.rds", scale))
}
out_rds <- file.path("inst", "extdata", "optop-benchmark-results.rds")
worker <- file.path("data-raw", "benchmark-worker.R")

# ---- corpora and fits (cached) ---------------------------------------------

simulate_corpus <- function(J, W, seed) {
  set.seed(seed)
  K_true <- 8L
  phi_true <- matrix(rgamma(K_true * W, 0.1), K_true, W)
  phi_true <- phi_true / rowSums(phi_true)
  counts <- matrix(0L, J, W,
                   dimnames = list(sprintf("doc%04d", seq_len(J)),
                                   sprintf("w%05d", seq_len(W))))
  for (j in seq_len(J)) {
    th <- rgamma(K_true, 0.4)
    th <- th / sum(th)
    p <- as.numeric(th %*% phi_true)
    counts[j, ] <- as.integer(rmultinom(1, sample(600:1500, 1), p))
  }
  empty <- which(colSums(counts) == 0)
  if (length(empty)) counts[1, empty] <- 1L
  counts
}

for (nm in names(scales)) {
  if (file.exists(fits_rds_for(nm))) {
    message("using cached fits in ", fits_rds_for(nm))
    next
  }
  sc <- scales[[nm]]
  message("fitting scale ", nm, " (J = ", sc$J, ", W = ", sc$W, ")")
  counts <- simulate_corpus(sc$J, sc$W, seed = 20260705 + sc$J)
  models <- lapply(K_grid, function(k) {
    message("  k = ", k)
    topicmodels::LDA(counts, k = k, method = "VEM",
                     control = list(seed = 500 + k,
                                    em = list(iter.max = em_iter)))
  })
  saveRDS(list(models = models,
               wdfm = quanteda::dfm_weight(quanteda::as.dfm(counts),
                                           scheme = "prop"),
               doc_lengths = rowSums(counts)),
          fits_rds_for(nm))
}

# ---- benchmark grid ---------------------------------------------------------

configs <- do.call(rbind, lapply(names(scales), function(nm) {
  rbind(
    data.frame(scale = nm, mode = "boot_r", threads = 1L),
    expand.grid(scale = nm, mode = c("stat", "boot"),
                threads = threads_sweep, stringsAsFactors = FALSE)
  )
}))

run_config <- function(scale, mode, threads) {
  res_file <- tempfile(fileext = ".rds")
  time_file <- tempfile(fileext = ".txt")
  status <- system2("/usr/bin/time",
                    args = c("-v", "-o", time_file, "Rscript", worker,
                             fits_rds_for(scale), scale, mode, threads,
                             res_file),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0L) stop("worker failed for ", scale, "/", mode, "/", threads)
  out <- readRDS(res_file)
  rss_line <- grep("Maximum resident set size", readLines(time_file),
                   value = TRUE)
  out$max_rss_mb <- as.numeric(sub(".*: ", "", rss_line)) / 1024
  out
}

runs <- vector("list", nrow(configs))
for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, ]
  message(sprintf("[%d/%d] %s / %s / %d thread(s)",
                  i, nrow(configs), cfg$scale, cfg$mode, cfg$threads))
  runs[[i]] <- run_config(cfg$scale, cfg$mode, cfg$threads)
}

results <- do.call(rbind, lapply(runs, function(r) {
  sc <- scales[[r$scale]]
  data.frame(scale = r$scale, J = sc$J, W = sc$W,
             n_models = length(K_grid), mode = r$mode, threads = r$threads,
             elapsed_sec = r$elapsed, max_rss_mb = r$max_rss_mb)
}))

# ---- correctness checks -----------------------------------------------------

pick <- function(scale, mode, threads) {
  for (r in runs) {
    if (r$scale == scale && r$mode == mode && r$threads == threads) {
      return(r$result)
    }
  }
  NULL
}

for (nm in names(scales)) {
  # thread sweeps are bit-identical
  for (mode in c("stat", "boot")) {
    base <- pick(nm, mode, 1L)
    for (th in setdiff(threads_sweep, 1L)) {
      stopifnot(identical(base, pick(nm, mode, th)))
    }
  }
  # the C++ bootstrap agrees with the R baseline within MC resolution:
  # different RNG streams, so p-values differ by O(1/sqrt(B))
  p_cpp <- pick(nm, "boot", 1L)$pval
  p_r <- pick(nm, "boot_r", 1L)$pval
  stopifnot(max(abs(p_cpp - p_r)) <= 5 * sqrt(0.25 / n_boot) + 2 / (n_boot + 1))
}
message("correctness checks passed")

# ---- store ------------------------------------------------------------------

meta <- list(
  generated = Sys.time(),
  cpu = tryCatch(
    sub("model name\\s*:\\s*", "",
        grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]),
    error = function(e) NA_character_
  ),
  cores = parallel::detectCores(),
  r_version = R.version.string,
  openmp = OpTop:::optop_openmp_available(),
  blas = extSoftVersion()[["BLAS"]],
  K_grid = K_grid, em_iter = em_iter, n_boot = n_boot, q = 0.95,
  scales = scales
)

saveRDS(list(results = results, meta = meta), out_rds)
message("written ", out_rds)
print(results)
