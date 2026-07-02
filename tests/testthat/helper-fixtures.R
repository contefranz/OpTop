# Shared fixtures for the OpTop test suite.
#
# A small deterministic corpus and a grid of LDA fits, built once per test run
# and cached in an environment so that every test file works on identical
# objects. Documents are long enough (relative to the vocabulary size) that the
# harmonized partition contains a genuine mix of rare and non-rare words, which
# keeps the index tests meaningful.

.optop_test_cache <- new.env(parent = emptyenv())

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
