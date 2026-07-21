# Shared fixtures for the OpTop test suite.
#
# A small deterministic corpus and a grid of LDA fits, built once per test run
# and cached in an environment so that every test file works on identical
# objects. Documents are long enough (relative to the vocabulary size) that the
# harmonized partition contains a genuine mix of rare and non-rare words, which
# keeps the index tests meaningful.

.optop_test_cache <- new.env(parent = emptyenv())

# Row-wise word proportions as a weighted quanteda dfm and as the plain dense
# matrix the reference implementations consume.
optop_wprop_fixture <- function(fx) {
  wdfm <- quanteda::dfm_weight(quanteda::as.dfm(fx$counts), scheme = "prop")
  W_prop <- sweep(fx$counts, 1, rowSums(fx$counts), "/")
  list(wdfm = wdfm, W_prop = W_prop)
}

# Gibbs fits on the same corpus, for the engine-generalization tests: same
# grid of K as the VEM fixture, cached alongside it.
optop_gibbs_fixture <- function(fx) {
  if (is.null(.optop_test_cache$gibbs)) {
    .optop_test_cache$gibbs <- lapply(2:4, function(k) {
      topicmodels::LDA(fx$counts, k = k, method = "Gibbs",
                       control = list(seed = 2000 + k, burnin = 100,
                                      iter = 300))
    })
  }
  .optop_test_cache$gibbs
}

# Iteration-capped CTM fits (slow to converge, but the adapters and the
# statistic only need valid fitted probabilities), cached alongside.
optop_ctm_fixture <- function(fx) {
  if (is.null(.optop_test_cache$ctm)) {
    .optop_test_cache$ctm <- lapply(2:3, function(k) {
      topicmodels::CTM(fx$counts, k = k,
                       control = list(seed = 3000 + k,
                                      em = list(iter.max = 30),
                                      var = list(iter.max = 20)))
    })
  }
  .optop_test_cache$ctm
}

optop_test_fixture <- function() {
  if (!is.null(.optop_test_cache$fixture)) {
    return(.optop_test_cache$fixture)
  }

  set.seed(20260702)
  J <- 30L   # documents
  W <- 120L  # vocabulary

  vocab <- sprintf("w%03d", seq_len(W))
  docs <- sprintf("doc%02d", seq_len(J))

  # Two latent word profiles mixed per document, so LDA has structure to find.
  prof1 <- rgamma(W, shape = ifelse(seq_len(W) <= W / 2, 4, 0.4))
  prof2 <- rgamma(W, shape = ifelse(seq_len(W) <= W / 2, 0.4, 4))
  prof1 <- prof1 / sum(prof1)
  prof2 <- prof2 / sum(prof2)

  counts <- matrix(0L, nrow = J, ncol = W, dimnames = list(docs, vocab))
  for (j in seq_len(J)) {
    mix <- rbeta(1, 2, 2)
    p <- mix * prof1 + (1 - mix) * prof2
    L_j <- sample(400:800, 1)
    counts[j, ] <- as.integer(rmultinom(1, size = L_j, prob = p))
  }
  # Guarantee no empty column (topicmodels needs every term observed).
  empty <- which(colSums(counts) == 0)
  if (length(empty)) counts[1, empty] <- counts[1, empty] + 1L

  models <- lapply(2:4, function(k) {
    topicmodels::LDA(counts, k = k, method = "VEM",
                     control = list(seed = 1000 + k))
  })

  dtm <- methods::as(counts, "CsparseMatrix")

  partition <- optop_make_partition(models, dtm, c = 5)
  baseline <- optop_make_baseline(dtm)

  fixture <- list(counts = counts, dtm = dtm, models = models,
                  partition = partition, baseline = baseline,
                  J = J, W = W)
  .optop_test_cache$fixture <- fixture
  fixture
}

# One model's Test 1 statistic and per-document envelope bins, through the
# compiled core — shared by the calibration and parallelization tests.
optop_envelope <- function(model, wdfm, q = 0.95) {
  doc_map <- seq_len(quanteda::ndoc(wdfm)) - 1L
  tp <- optop_as_theta_phi(model)
  dfm_t <- Matrix::t(methods::as(wdfm, "dgCMatrix"))
  out <- optimal_topic_core(tp$theta, tp$phi, dfm_t@p, dfm_t@i, dfm_t@x,
                            nrow(dfm_t), q, doc_map,
                            return_envelope = TRUE, n_threads = 1L)
  list(stat = out$stat,
       probs = .optop_split_envelope(out$bin_probs, out$bin_counts))
}

# Train/eval split fixture for the held-out and moment-test suites: models
# fitted on the first 20 documents of the shared corpus, evaluated on the
# remaining 10, with the TRAINING baseline.
optop_holdout_fixture <- function() {
  if (is.null(.optop_test_cache$holdout)) {
    fx <- optop_test_fixture()
    tr <- fx$dtm[1:20, ]
    ev <- fx$dtm[21:30, ]
    models <- lapply(2:4, function(k) {
      topicmodels::LDA(as.matrix(tr), k = k, method = "VEM",
                       control = list(seed = 100 + k))
    })
    .optop_test_cache$holdout <- list(tr = tr, ev = ev, models = models,
                                      baseline = optop_make_baseline(tr))
  }
  .optop_test_cache$holdout
}
