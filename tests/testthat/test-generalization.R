# Generalizability (v0.11.0): the optop_as_theta_phi() adapter layer, the
# model-agnostic C++ core, and the identifier-based alignment contract.
# topicmodels engines (VEM, Gibbs, CTM) are characterized against the naive
# Eq. (8) reference through their gamma/beta slots; seededlda and the
# NLPstudio-style wrapper against the same reference on plain (theta, phi)
# pairs.

run_ot <- function(...) {
  suppressMessages(optimal_topic(..., do_plot = FALSE))
}

# an NLPstudio-shaped fit built by hand from any topicmodels model, so the
# nlp_topic_fit adapter is exercised without depending on NLPstudio
nlp_wrap <- function(model) {
  p <- topicmodels::posterior(model)
  structure(list(engine = "topicmodels", model = "lda", method = "VEM",
                 model_object = NULL,
                 dtw = p$topics, tww = p$terms,
                 doc_ids = as.character(model@documents),
                 vocab = as.character(model@terms)),
            class = c("nlp_topic_fit", "list"))
}

test_that("adapters return the full contract for every supported engine", {
  fx <- optop_test_fixture()

  fits <- c(fx$models[1L],                      # LDA_VEM
            optop_gibbs_fixture(fx)[1L],        # LDA_Gibbs
            optop_ctm_fixture(fx)[1L],          # CTM_VEM
            list(nlp_wrap(fx$models[[1L]])))    # nlp_topic_fit
  for (fit in fits) {
    tp <- optop_as_theta_phi(fit)
    expect_named(tp, c("theta", "phi", "K", "docs", "terms"))
    expect_identical(nrow(tp$theta), fx$J)
    expect_identical(ncol(tp$phi), fx$W)
    expect_identical(ncol(tp$theta), nrow(tp$phi))
    expect_identical(tp$docs, rownames(fx$counts))
    expect_identical(tp$terms, colnames(fx$counts))
    expect_true(.optop_validate_theta_phi(tp, class(fit)))
  }
})

test_that("unsupported classes fail with the supported-engine list", {
  expect_error(optop_as_theta_phi(lm(dist ~ speed, cars)),
               "unsupported topic model class <lm>")
})

test_that("a Gibbs grid matches the naive Eq. (8) reference", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  gibbs <- optop_gibbs_fixture(fx)

  got <- run_ot(gibbs, wp$wdfm)
  ref <- ref_optimal_topic(gibbs, wp$W_prop, q = 0.95)
  expect_equal(got$topic, ref$topic)
  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10)
  expect_equal(got$df, ref$df)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10)
})

test_that("tied fitted probabilities resolve like R's stable order()", {
  # exact ties by construction: every vocabulary column duplicated, so the
  # envelope boundary straddles tied values in essentially every document.
  # The compiled comparator (probability descending, index ascending) must
  # reproduce order(X, decreasing = TRUE) exactly on every platform.
  set.seed(42)
  K <- 3L; W <- 40L; J <- 12L
  phi_half <- matrix(rgamma(K * (W / 2), 1), K, W / 2)
  phi <- cbind(phi_half, phi_half)
  phi <- phi / rowSums(phi)
  theta <- matrix(rgamma(J * K, 1), J, K)
  theta <- theta / rowSums(theta)
  colnames(phi) <- sprintf("w%02d", seq_len(W))
  rownames(theta) <- sprintf("d%02d", seq_len(J))
  o <- matrix(rgamma(J * W, 1), J, W)
  o <- o / rowSums(o)

  q <- 0.95
  dfm_t <- Matrix::t(methods::as(Matrix::Matrix(o, sparse = TRUE),
                                 "CsparseMatrix"))
  got <- optimal_topic_core(theta, phi, dfm_t@p, dfm_t@i, dfm_t@x,
                            nrow(dfm_t), q,
                            doc_map = seq_len(J) - 1L,
                            return_envelope = TRUE, n_threads = 1L)

  X <- theta %*% phi
  chi <- df <- numeric(J)
  bins_ref <- vector("list", J)
  for (j in seq_len(J)) {
    ord <- order(X[j, ], decreasing = TRUE)   # stable: ties by index
    x <- X[j, ord]
    oo <- o[j, ord]
    p_j <- which(cumsum(x) > q)[1]
    if (is.na(p_j)) p_j <- W
    head_idx <- seq_len(p_j)
    pearson <- sum((oo[head_idx] - x[head_idx])^2 / x[head_idx])
    n_bins <- p_j
    bins_ref[[j]] <- x[head_idx]
    if (p_j < W) {
      tail_x <- sum(x) - sum(x[head_idx])
      tail_o <- sum(oo) - sum(oo[head_idx])
      pearson <- pearson + (tail_o - tail_x)^2 / tail_x
      n_bins <- n_bins + 1L
      bins_ref[[j]] <- c(bins_ref[[j]], tail_x)
    }
    chi[j] <- n_bins * pearson
    df[j] <- n_bins - 1L
  }

  expect_equal(got$stat[1, 2], sum(chi), tolerance = 1e-12)
  expect_equal(got$stat[1, 3], sum(df))
  expect_equal(OpTop:::.optop_split_envelope(got$bin_probs, got$bin_counts),
               bins_ref, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("a CTM grid matches the naive Eq. (8) reference", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  ctm <- optop_ctm_fixture(fx)

  got <- run_ot(ctm, wp$wdfm)
  ref <- ref_optimal_topic(ctm, wp$W_prop, q = 0.95)
  expect_equal(got$topic, ref$topic)
  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10)
})

test_that("engines can be mixed within one grid", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  gibbs <- optop_gibbs_fixture(fx)

  mixed <- list(fx$models[[1L]], gibbs[[2L]], fx$models[[3L]])  # k = 2, 3, 4
  got <- run_ot(mixed, wp$wdfm)
  ref <- ref_optimal_topic(mixed, wp$W_prop, q = 0.95)
  expect_equal(got$topic, c(2, 3, 4))
  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10)
})

test_that("all three seededlda constructors run and match the reference", {
  skip_if_not_installed("seededlda")
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  dfmat <- quanteda::as.dfm(fx$counts)

  dict <- quanteda::dictionary(list(low = c("w001", "w002", "w003"),
                                    high = c("w061", "w062", "w063")))
  set.seed(11)
  seeded <- seededlda::textmodel_seededlda(dfmat, dict, max_iter = 200,
                                           verbose = FALSE)
  plain <- seededlda::textmodel_lda(dfmat, k = 3, max_iter = 200,
                                    verbose = FALSE)
  # seqlda warns that gamma is idle when documents are not split into
  # sentences; irrelevant here, the fit is what matters
  seq <- suppressWarnings(
    seededlda::textmodel_seqlda(dfmat, k = 4, max_iter = 200, verbose = FALSE)
  )

  grid <- list(seeded, plain, seq)  # k = 2, 3, 4, one per constructor
  got <- run_ot(grid, wp$wdfm)
  tps <- lapply(grid, function(f) list(theta = f$theta, phi = f$phi))
  ref <- ref_optimal_topic_tp(tps, wp$W_prop, q = 0.95)
  expect_equal(got$topic, c(2, 3, 4))
  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10)
  expect_equal(got$df, ref$df)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10)
})

test_that("an nlp_topic_fit grid reproduces the raw-model results exactly", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  wrapped <- lapply(fx$models, nlp_wrap)
  got <- run_ot(wrapped, wp$wdfm)
  ref <- run_ot(fx$models, wp$wdfm)
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("nlp_topic_fit falls back to model_object and errors when empty", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  hollow <- lapply(fx$models, function(m) {
    structure(list(dtw = NULL, tww = NULL, model_object = m),
              class = c("nlp_topic_fit", "list"))
  })
  got <- run_ot(hollow, wp$wdfm)
  ref <- run_ot(fx$models, wp$wdfm)
  expect_equal(got, ref, tolerance = 1e-12)

  empty <- structure(list(dtw = NULL, tww = NULL, model_object = NULL),
                     class = c("nlp_topic_fit", "list"))
  expect_error(optop_as_theta_phi(empty), "stores neither")
})

test_that("alignment is per model: document order inside a model is irrelevant", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  # one model of the grid saw the documents in a different order
  set.seed(7)
  perm <- sample(fx$J)
  shuffled <- nlp_wrap(fx$models[[2L]])
  shuffled$dtw <- shuffled$dtw[perm, , drop = FALSE]
  shuffled$doc_ids <- shuffled$doc_ids[perm]

  grid <- list(nlp_wrap(fx$models[[1L]]), shuffled, nlp_wrap(fx$models[[3L]]))
  got <- run_ot(grid, wp$wdfm)
  ref <- run_ot(fx$models, wp$wdfm)
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("documents missing from any model are dropped everywhere", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  # the second model never saw the first document
  partial <- nlp_wrap(fx$models[[2L]])
  partial$dtw <- partial$dtw[-1L, , drop = FALSE]
  partial$doc_ids <- partial$doc_ids[-1L]
  grid <- list(nlp_wrap(fx$models[[1L]]), partial, nlp_wrap(fx$models[[3L]]))

  expect_message(
    optimal_topic(grid, wp$wdfm, do_plot = FALSE, verbose = FALSE),
    "Removed 1 document"
  )

  # equivalent to evaluating the full grid on the reduced dfm
  got <- run_ot(grid, wp$wdfm)
  ref <- run_ot(fx$models, wp$wdfm[-1L, ])
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("models without document identifiers are rejected", {
  fx <- optop_test_fixture()

  anon <- nlp_wrap(fx$models[[1L]])
  anon$doc_ids <- NULL
  rownames(anon$dtw) <- NULL
  expect_error(.optop_validate_theta_phi(optop_as_theta_phi(anon),
                                         class(anon)),
               "must expose one document identifier")

  dup <- nlp_wrap(fx$models[[1L]])
  dup$doc_ids[2L] <- dup$doc_ids[1L]
  expect_error(.optop_validate_theta_phi(optop_as_theta_phi(dup),
                                         class(dup)),
               "duplicated document identifiers")
})

test_that("feature order is realigned, other mismatches are errors", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  # same features, different order: fixed with a signal, results invariant
  set.seed(9)
  perm <- sample(fx$W)
  wdfm_perm <- wp$wdfm[, perm]
  expect_message(
    optimal_topic(fx$models, wdfm_perm, do_plot = FALSE, verbose = FALSE),
    "Reordered the features"
  )
  expect_equal(run_ot(fx$models, wdfm_perm),
               run_ot(fx$models, wp$wdfm), tolerance = 1e-12)

  # a missing feature cannot be fixed
  expect_error(run_ot(fx$models, wp$wdfm[, -1L]),
               "do not match the models' vocabulary")

  # models disagreeing on the vocabulary cannot be compared
  alien <- nlp_wrap(fx$models[[2L]])
  alien$vocab[1L] <- "not-a-word"
  colnames(alien$tww)[1L] <- "not-a-word"
  expect_error(run_ot(list(nlp_wrap(fx$models[[1L]]), alien), wp$wdfm),
               "share the same vocabulary")
})

test_that("the grid is sorted by K and duplicated K is an error", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_message(
    optimal_topic(rev(fx$models), wp$wdfm, do_plot = FALSE, verbose = FALSE),
    "Reordered topic_models"
  )
  got <- run_ot(rev(fx$models), wp$wdfm)
  ref <- run_ot(fx$models, wp$wdfm)
  expect_equal(got, ref, tolerance = 1e-12)

  expect_error(run_ot(fx$models[c(1L, 1L, 2L)], wp$wdfm),
               "more than one model with the same number of topics")
})

test_that("a non-proportion dfm is flagged", {
  fx <- optop_test_fixture()

  # selection = "min" keeps the sequential fallback alert out of the way:
  # a counts dfm rejects everywhere, which is exactly the point of the flag
  expect_message(
    optimal_topic(fx$models, quanteda::as.dfm(fx$counts), selection = "min",
                  do_plot = FALSE, verbose = FALSE),
    "rows do not sum to 1"
  )
})

test_that("the deprecated lda_models argument still works and warns", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_warning(
    got <- run_ot(lda_models = fx$models, weighted_dfm = wp$wdfm),
    class = "lifecycle_warning_deprecated"
  )
  ref <- run_ot(topic_models = fx$models, weighted_dfm = wp$wdfm)
  expect_identical(got, ref)
})

test_that("calibration runs on non-VEM engines", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  gibbs <- optop_gibbs_fixture(fx)
  doc_lens <- stats::setNames(rowSums(fx$counts), rownames(fx$counts))

  got <- run_ot(gibbs, wp$wdfm, calibrate = "bootstrap", n_boot = 50,
                doc_lengths = doc_lens, seed = 99)
  expect_identical(names(got),
                   c("topic", "OpTop", "df", "pval", "pval_chisq"))
  expect_true(all(got$pval >= 0 & got$pval <= 1))

  got_mm <- run_ot(gibbs, wp$wdfm, calibrate = "moment",
                   doc_lengths = doc_lens)
  expect_true(all(got_mm$pval >= 0 & got_mm$pval <= 1))
})
