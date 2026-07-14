# The family-generic null-discrepancy floor (v0.14.1): documents with
# D_j(null) < min_null leave the document-level aggregation, with default
# min_null = partition$c. The three-regime corpus below reproduces the
# pathology in miniature:
# * docs A, B: long and far from the pooled baseline (healthy);
# * doc C: L = 40 tokens over W = 400 words, so tau_j = c/L > every fitted
#   and baseline probability: the whole support is swept into the min bin,
#   observed == baseline by construction, D_null ~ 0 in every family at
#   once (the structural-collapse regime, exact 0/0);
# * doc D: one token away from the pooled corpus distribution, so
#   0 < D_null < c (the near-baseline instability regime).

null_floor_fixture <- function() {
  set.seed(11)
  W <- 400
  vocab <- sprintf("w%03d", seq_len(W))
  A <- stats::rmultinom(1, 3000, prob = c(rep(4, 100), rep(1, 300)))[, 1]
  B <- stats::rmultinom(1, 3000, prob = c(rep(1, 300), rep(4, 100)))[, 1]
  C <- stats::rmultinom(1, 40, prob = rep(1, W))[, 1]
  D <- A + B + C
  D[1] <- D[1] + 1
  D[2] <- D[2] - 1
  N <- rbind(A = A, B = B, C = C, D = D)
  colnames(N) <- vocab
  dtm <- methods::as(Matrix::Matrix(N, sparse = TRUE), "CsparseMatrix")
  theta <- matrix(1, 4, 1, dimnames = list(rownames(N), NULL))
  phi <- matrix(1 / W, 1, W, dimnames = list(NULL, vocab))
  models <- list(OpTop:::.optop_tp(theta, phi))
  list(dtm = dtm, models = models,
       partition = optop_make_partition(models, dtm, c = 1),
       baseline = optop_make_baseline(dtm))
}

index_fun <- function(metric) {
  switch(metric, se = optop_index_se, chisq = optop_index_chisq,
         deviance = optop_index_deviance)
}

test_that("the default floor excludes both degenerate regimes, all families", {
  fx <- null_floor_fixture()
  for (metric in c("se", "chisq", "deviance")) {
    res <- suppressMessages(
      index_fun(metric)(fx$models[[1]], fx$dtm, fx$partition, fx$baseline,
                        macro = TRUE)
    )
    # doc C collapses structurally, doc D sits below the floor
    expect_lt(abs(res$d_null[3]), 1e-6)
    expect_gt(res$d_null[4], 0)
    expect_lt(res$d_null[4], 1)
    expect_identical(unname(is.na(res$r2_doc)), c(FALSE, FALSE, TRUE, TRUE))
    expect_identical(res$n_null_excluded, 2L)
    expect_equal(res$null_excluded_share, 0.5)
    # Macro over the kept documents only; Micro sums restricted to J+
    expect_equal(res$r2_macro, mean(res$r2_doc[1:2]))
    expect_equal(res$r2,
                 1 - sum(res$d_model[1:2]) / sum(res$d_null[1:2]))
    # independent reference computation with the same floor
    ref <- ref_index_document(fx$models[[1]], fx$dtm, fx$partition,
                              fx$baseline, metric = metric, min_null = 1)
    expect_equal(res$r2, ref$r2, tolerance = 1e-10)
    expect_equal(res$r2_macro, ref$r2_macro, tolerance = 1e-10)
    expect_identical(res$n_null_excluded, ref$n_null_excluded)
  }
})

test_that("min_null = 0 reproduces the legacy strict-positivity rule", {
  # on collapsed documents the legacy ratio amplifies last-ulp roundoff of
  # D_null ~ 1e-28 beyond any cross-implementation tolerance (the pathology
  # itself), so the legacy path is verified as a package-internal identity
  # on the returned discrepancies rather than against the slow reference
  fx <- null_floor_fixture()
  for (metric in c("se", "chisq", "deviance")) {
    res <- index_fun(metric)(fx$models[[1]], fx$dtm, fx$partition,
                             fx$baseline, macro = TRUE, min_null = 0)
    valid <- res$d_null > 0
    r_exp <- rep(NA_real_, length(res$d_null))
    r_exp[valid] <- 1 - res$d_model[valid] / res$d_null[valid]
    expect_identical(res$r2_doc, r_exp)
    expect_identical(res$r2,
                     1 - sum(res$d_model[valid]) / sum(res$d_null[valid]))
    expect_identical(res$r2_macro, mean(r_exp[valid]))
    # the sub-floor documents destroy the unfloored Macro
    expect_lt(res$r2_macro, -100)
  }
})

test_that("healthy corpora are untouched by the default floor", {
  hfx <- optop_test_fixture()
  part <- optop_make_partition(hfx$models, hfx$dtm, c = 1)
  base <- optop_make_baseline(hfx$dtm)
  for (metric in c("se", "chisq", "deviance")) {
    res1 <- suppressMessages(
      index_fun(metric)(hfx$models[[2]], hfx$dtm, part, base, macro = TRUE)
    )
    res0 <- suppressMessages(
      index_fun(metric)(hfx$models[[2]], hfx$dtm, part, base, macro = TRUE,
                        min_null = 0)
    )
    expect_identical(res1$r2, res0$r2)
    expect_identical(res1$r2_macro, res0$r2_macro)
    expect_identical(res1$r2_doc, res0$r2_doc)
    expect_identical(res1$n_null_excluded, 0L)
  }
})

test_that("the Micro index moves far less than the Macro index", {
  fx <- null_floor_fixture()
  for (metric in c("se", "chisq", "deviance")) {
    res1 <- suppressMessages(
      index_fun(metric)(fx$models[[1]], fx$dtm, fx$partition, fx$baseline,
                        macro = TRUE)
    )
    res0 <- index_fun(metric)(fx$models[[1]], fx$dtm, fx$partition,
                              fx$baseline, macro = TRUE, min_null = 0)
    d_micro <- abs(res1$r2 - res0$r2)
    d_macro <- abs(res1$r2_macro - res0$r2_macro)
    expect_lt(d_micro, d_macro / 100)
  }
})

test_that("exclusions are reported once per call", {
  fx <- null_floor_fixture()
  expect_message(
    optop_index_deviance(fx$models[[1]], fx$dtm, fx$partition, fx$baseline),
    "null-discrepancy floor"
  )
  # a healthy corpus stays silent
  hfx <- optop_test_fixture()
  part <- optop_make_partition(hfx$models, hfx$dtm, c = 1)
  base <- optop_make_baseline(hfx$dtm)
  expect_silent(
    optop_index_deviance(hfx$models[[2]], hfx$dtm, part, base)
  )
})

test_that("optop_index_table applies and reports the floor", {
  fx <- null_floor_fixture()
  expect_message(
    tab <- optop_index_table(fx$models, fx$dtm,
                             partition = fx$partition,
                             baseline = fx$baseline, macro = TRUE),
    "null-discrepancy floor"
  )
  ref <- ref_index_document(fx$models[[1]], fx$dtm, fx$partition,
                            fx$baseline, metric = "deviance", min_null = 1)
  expect_equal(tab$R2_dev, ref$r2, tolerance = 1e-10)
  expect_equal(tab$R2_dev_macro, ref$r2_macro, tolerance = 1e-10)
  # the escape hatch restores the legacy table
  tab0 <- optop_index_table(fx$models, fx$dtm,
                            partition = fx$partition,
                            baseline = fx$baseline, macro = TRUE,
                            min_null = 0)
  ref0 <- ref_index_document(fx$models[[1]], fx$dtm, fx$partition,
                             fx$baseline, metric = "deviance", min_null = 0)
  expect_equal(tab0$R2_dev_macro, ref0$r2_macro, tolerance = 1e-10)
})

test_that("min_null is validated", {
  fx <- null_floor_fixture()
  expect_error(
    optop_index_deviance(fx$models[[1]], fx$dtm, fx$partition, fx$baseline,
                         min_null = -1),
    "nonnegative"
  )
  expect_error(
    optop_index_deviance(fx$models[[1]], fx$dtm, fx$partition, fx$baseline,
                         min_null = c(1, 2)),
    "single"
  )
  expect_error(
    optop_index_holdout(fx$models, fx$dtm, fx$baseline, min_null = Inf),
    "nonnegative|single"
  )
})

# Held-out three-regime construction: models are fitted on a training split
# and the evaluation set plants the same two degenerate regimes against the
# TRAINING baseline (doc "dup" duplicates the pooled training counts, so its
# binned distribution equals pi_tr exactly; doc "tiny" collapses its support).
null_floor_holdout_fixture <- function() {
  if (is.null(.optop_test_cache$null_floor_ho)) {
    set.seed(23)
    W <- 120
    vocab <- sprintf("v%03d", seq_len(W))
    N_tr <- t(vapply(seq_len(12), function(j) {
      pr <- rep(1, W)
      pr[((j - 1) %% 4) * 30 + seq_len(30)] <- 6
      stats::rmultinom(1, 600, prob = pr)[, 1]
    }, numeric(W)))
    dimnames(N_tr) <- list(paste0("tr", seq_len(12)), vocab)
    models <- lapply(2:3, function(k) {
      topicmodels::LDA(N_tr, k = k, method = "VEM",
                       control = list(seed = 400 + k))
    })
    ev_healthy <- t(vapply(seq_len(4), function(j) {
      pr <- rep(1, W)
      pr[((j - 1) %% 4) * 30 + seq_len(30)] <- 6
      stats::rmultinom(1, 500, prob = pr)[, 1]
    }, numeric(W)))
    tiny <- stats::rmultinom(1, 12, prob = rep(1, W))[, 1]
    # one token away from the pooled training counts: the binned distance to
    # pi_tr is strictly positive yet far below the floor for the Pearson and
    # deviance families (the squared-error distance is exactly 2)
    dup <- colSums(N_tr)
    dup[1] <- dup[1] + 1
    dup[2] <- dup[2] - 1
    N_ev <- rbind(ev_healthy, tiny = tiny, dup = dup)
    rownames(N_ev) <- c(paste0("ev", seq_len(4)), "tiny", "dup")
    dtm_eval <- methods::as(Matrix::Matrix(N_ev, sparse = TRUE),
                            "CsparseMatrix")
    .optop_test_cache$null_floor_ho <- list(
      models = models, dtm_eval = dtm_eval,
      baseline = optop_make_baseline(
        methods::as(Matrix::Matrix(N_tr, sparse = TRUE), "CsparseMatrix")
      )
    )
  }
  .optop_test_cache$null_floor_ho
}

test_that("the holdout floor matches the reference and repairs Macro", {
  fx <- null_floor_holdout_fixture()
  ho <- suppressMessages(
    optop_index_holdout(fx$models, fx$dtm_eval, fx$baseline)
  )
  ref <- ref_holdout(fx$models, fx$dtm_eval, fx$baseline, min_null = 1)
  for (metric in ho$metrics) {
    sm <- as.data.frame(ho$summary)
    sm <- sm[sm$metric == metric, ]
    for (i in seq_along(fx$models)) {
      r <- ref[[metric]][[i]]
      expect_gte(r$n_excluded, 1)
      expect_equal(sm$macro[i], r$macro, tolerance = 1e-8)
      expect_equal(sm$macro_se[i], r$macro_se, tolerance = 1e-8)
      expect_equal(sm$micro[i], r$micro, tolerance = 1e-8)
      expect_equal(sm$gap_se[i], r$gap_se, tolerance = 1e-8)
      expect_identical(sm$n_docs[i], r$n)
      expect_identical(sm$n_null_excluded[i], r$n_excluded)
      expect_identical(unname(is.na(ho$scores[[metric]][, i])),
                       unname(is.na(r$r)))
    }
  }
  # legacy escape hatch: package-internal identity on the returned
  # discrepancies (cross-implementation conformance is meaningless on
  # collapsed documents, where the ratio amplifies last-ulp roundoff)
  ho0 <- optop_index_holdout(fx$models, fx$dtm_eval, fx$baseline,
                             min_null = 0)
  for (metric in ho0$metrics) {
    sm0 <- as.data.frame(ho0$summary)
    sm0 <- sm0[sm0$metric == metric, ]
    for (i in seq_along(fx$models)) {
      D_K <- ho0$d_model[[metric]][, i]
      D_0 <- ho0$d_null[[metric]][, i]
      valid <- D_0 > 0
      expect_identical(sm0$macro[i],
                       mean(1 - D_K[valid] / D_0[valid]))
      expect_identical(sm0$n_docs[i], sum(valid))
    }
  }
  # the legacy rule keeps the sub-floor near-baseline document, whose score
  # 1 - D_K/D_null with D_null << 1 is unbounded in either direction, so
  # the Macro estimand shifts materially relative to the floored one
  sm0 <- as.data.frame(ho0$summary)
  sm1 <- as.data.frame(ho$summary)
  dev0 <- sm0[sm0$metric == "deviance", ]
  dev1 <- sm1[sm1$metric == "deviance", ]
  expect_true(all(dev0$n_docs > dev1$n_docs))
  expect_gt(max(abs(dev0$macro - dev1$macro)), 0.1)
})

test_that("the gain table drops floored documents pairwise", {
  fx <- null_floor_holdout_fixture()
  ho <- suppressMessages(
    optop_index_holdout(fx$models, fx$dtm_eval, fx$baseline)
  )
  gt <- optop_gain_table(ho, metric = "deviance")
  n_kept <- ho$summary$n_docs[ho$summary$metric == "deviance"][1]
  expect_identical(gt$gains$n, n_kept)
})

test_that("the floor composes with stabilization", {
  fx <- null_floor_holdout_fixture()
  # min_null = 0 with stabilization keeps every document (legacy semantics)
  ho_leg <- optop_index_holdout(fx$models, fx$dtm_eval, fx$baseline,
                                min_null = 0, stabilize = 1)
  expect_true(all(ho_leg$summary$n_docs == nrow(fx$dtm_eval)))
  expect_true(all(is.finite(ho_leg$scores$deviance)))
  # with the floor, stabilization transforms only the kept denominators
  ho_st <- suppressMessages(
    optop_index_holdout(fx$models, fx$dtm_eval, fx$baseline,
                        min_null = 1, stabilize = 5)
  )
  D_K <- ho_st$d_model$deviance[, 1]
  D_0 <- ho_st$d_null$deviance[, 1]
  keep <- D_0 >= 1
  manual <- rep(NA_real_, length(D_K))
  manual[keep] <- 1 - D_K[keep] / pmax(D_0[keep], 5)
  expect_equal(unname(ho_st$scores$deviance[, 1]), unname(manual),
               tolerance = 1e-12)
})
