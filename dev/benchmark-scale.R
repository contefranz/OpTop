# Scale benchmark for the 0.15.0 architecture: sparse partition, merge-join
# index kernels, zero-copy boundary, and sharded evaluation. Records
# reference timings on synthetic Zipf corpora at sizes the pre-0.15.0
# dense-mask code path could not represent (a dense J x W rare mask at the
# largest configuration below would need 4 * J * W bytes = 200 GB).
#
# Run from the package root: Rscript dev/benchmark-scale.R
# Output: dev/benchmark-scale.md

library(OpTop)
library(Matrix)

n_threads <- max(1L, parallel::detectCores() - 2L)

zipf_corpus <- function(J, W, distinct, mean_count, seed) {
  set.seed(seed)
  zipf <- 1 / seq_len(W)
  jj <- sample.int(W, J * distinct, replace = TRUE, prob = zipf)
  ii <- rep.int(seq_len(J), rep(distinct, J))
  x <- stats::rpois(length(jj), mean_count) + 1
  m <- Matrix::sparseMatrix(i = ii, j = jj, x = x, dims = c(J, W),
                            dimnames = list(sprintf("d%07d", seq_len(J)),
                                            sprintf("w%05d", seq_len(W))))
  methods::as(m, "CsparseMatrix")
}

zdir <- function(k, W, zipf_w) {
  g <- matrix(stats::rgamma(k * W, shape = 0.3), k) *
    matrix(zipf_w, k, W, byrow = TRUE)
  g / rowSums(g)
}

rdir <- function(n, k) {
  g <- matrix(stats::rgamma(n * k, 1), n)
  g / rowSums(g)
}

toy_models <- function(dtm, ks, seed) {
  set.seed(seed)
  J <- nrow(dtm); W <- ncol(dtm)
  zipf_w <- 1 / seq_len(W)
  lapply(ks, function(k) {
    theta <- rdir(J, k)
    rownames(theta) <- rownames(dtm)
    phi <- zdir(k, W, zipf_w)
    colnames(phi) <- colnames(dtm)
    OpTop:::.optop_tp(theta, phi)
  })
}

bench_one <- function(J, W, distinct, mean_count, ks, shards, seed) {
  dtm <- zipf_corpus(J, W, distinct, mean_count, seed)
  models <- toy_models(dtm, ks, seed + 1)
  idx <- parallel::splitIndices(J, shards)
  corp <- optop_corpus(lapply(idx, function(i) dtm[i, , drop = FALSE]))

  t_part <- system.time(
    part <- optop_make_partition(models, corp, c = 1, n_threads = n_threads)
  )[3]
  base <- optop_make_baseline(corp)
  t_doc <- system.time(
    res <- suppressMessages(
      optop_index_deviance(models[[length(models)]], corp, part, base,
                           macro = TRUE, n_threads = n_threads)
    )
  )[3]
  # row-normalized proportions per shard for the Test 1 sweep (the diagonal
  # product drops dimnames, so they are restored per shard)
  wdfm_sh <- lapply(idx, function(i) {
    m <- dtm[i, , drop = FALSE]
    w <- methods::as(Matrix::Diagonal(x = 1 / Matrix::rowSums(m)) %*% m,
                     "CsparseMatrix")
    dimnames(w) <- dimnames(m)
    w
  })
  wcorp <- optop_corpus(wdfm_sh)
  t_ot <- system.time(
    ot <- optimal_topic(models, wcorp, do_plot = FALSE, verbose = FALSE,
                        n_threads = n_threads)
  )[3]

  data.frame(J = J, W = W, nnz = length(dtm@x), shards = shards,
             ks = paste(ks, collapse = "/"),
             partition_s = round(unname(t_part), 1),
             partition_mb = round(as.numeric(object.size(part)) / 1024^2, 1),
             doc_index_s = round(unname(t_doc), 1),
             optimal_topic_s = round(unname(t_ot), 1),
             micro_dev = round(res$r2, 4),
             k_pick = ot$topic[which.min(ot$OpTop)])
}

cat("threads:", n_threads, "\n")
rows <- rbind(
  bench_one(100000L, 30000L, 60L, 30, c(5L, 8L), shards = 1L, seed = 1),
  bench_one(500000L, 50000L, 60L, 30, c(5L, 8L), shards = 4L, seed = 2),
  bench_one(1000000L, 50000L, 50L, 30, c(5L, 8L), shards = 4L, seed = 3)
)
print(rows, row.names = FALSE)

md <- c(
  "# Scale benchmark (0.15.0 architecture)",
  "",
  sprintf("Threads: %d. Synthetic Zipf corpora; two toy models per grid.",
          n_threads),
  "A dense J x W rare mask at the largest size would need 200 GB;",
  "the sparse partition sizes below are the complete objects.",
  "",
  knitr::kable(rows, format = "markdown")
)
writeLines(md, file.path("dev", "benchmark-scale.md"))
cat("written dev/benchmark-scale.md\n")
