# The held-out layer of the paper's Section 3.7: the package chain
# (alignment, OOV convention, fold-in, held-out partition with the training
# baseline, discrepancies, Macro/Micro/gap inference) must reproduce the
# explicit-loop reference, and the epsilon-adequacy rule must follow
# Definition 1 exactly.

test_that("the held-out chain matches the naive reference", {
  hx <- optop_holdout_fixture()
  ho <- optop_index_holdout(hx$models, hx$ev, hx$baseline, c = 1)
  ref <- ref_holdout(hx$models, hx$ev, hx$baseline, c = 1)

  for (metric in c("se", "chisq", "deviance")) {
    for (i in seq_along(hx$models)) {
      r <- ref[[metric]][[i]]
      sm <- as.data.frame(ho$summary)
      row <- sm[sm$metric == metric, ][i, ]
      info <- sprintf("metric=%s K=%d", metric, ho$K[i])
      expect_equal(unname(ho$scores[[metric]][, i]), unname(r$r),
                   tolerance = 1e-10, info = info)
      expect_equal(row$macro, r$macro, tolerance = 1e-10, info = info)
      expect_equal(row$macro_se, r$macro_se, tolerance = 1e-10, info = info)
      expect_equal(c(row$ci_lo, row$ci_hi), r$ci,
                   tolerance = 1e-10, info = info)
      expect_equal(row$micro, r$micro, tolerance = 1e-10, info = info)
      expect_equal(row$micro_se, r$micro_se, tolerance = 1e-10, info = info)
      expect_equal(row$gap, r$gap, tolerance = 1e-10, info = info)
      expect_equal(row$gap_se, r$gap_se, tolerance = 1e-10, info = info)
      expect_identical(row$n_docs, r$n)
      expect_equal(unname(ho$d_model[[metric]][, i]), unname(r$D_K),
                   tolerance = 1e-10, info = info)
      expect_equal(unname(ho$d_null[[metric]][, i]), unname(r$D_0),
                   tolerance = 1e-10, info = info)
    }
  }
})

test_that("out-of-support words follow the min-bin convention", {
  hx <- optop_holdout_fixture()
  # plant an evaluation-only word with visible mass
  extra <- Matrix::Matrix(c(5, rep(0, 9)), ncol = 1, sparse = TRUE,
                          dimnames = list(rownames(hx$ev), "zzz_new_word"))
  ev_oov <- methods::as(cbind(hx$ev, extra), "CsparseMatrix")

  expect_message(
    ho <- optop_index_holdout(hx$models, ev_oov, hx$baseline, c = 1),
    "out-of-support"
  )
  ref <- ref_holdout(hx$models, ev_oov, hx$baseline, c = 1)
  for (metric in c("se", "chisq", "deviance")) {
    expect_equal(unname(ho$scores[[metric]][, 2]),
                 unname(ref[[metric]][[2]]$r), tolerance = 1e-10)
  }

  # the planted token changes the affected document only
  ho0 <- optop_index_holdout(hx$models, hx$ev, hx$baseline, c = 1)
  d_new <- ho$d_null$deviance[, 1]
  d_old <- ho0$d_null$deviance[, 1]
  expect_false(isTRUE(all.equal(d_new[1], d_old[1])))
  expect_equal(unname(d_new[-1]), unname(d_old[-1]), tolerance = 1e-12)
})

test_that("evaluation documents without aligned tokens are dropped", {
  hx <- optop_holdout_fixture()
  ev <- as.matrix(hx$ev)
  ev[3, ] <- 0
  expect_warning(
    ho <- optop_index_holdout(hx$models, methods::as(ev, "CsparseMatrix"),
                              hx$baseline, c = 1, metrics = "deviance"),
    "no token"
  )
  expect_identical(nrow(ho$scores$deviance), nrow(ev) - 1L)
})

test_that("the stabilized index makes every document score finite", {
  hx <- optop_holdout_fixture()
  ho <- optop_index_holdout(hx$models, hx$ev, hx$baseline, c = 1,
                            metrics = "deviance", stabilize = 1e-8)
  expect_false(anyNA(ho$scores$deviance))
  # r_j = 1 - D_K / max(D_null, epsilon), here D_null > epsilon everywhere,
  # so the stabilized scores coincide with the plain ones
  plain <- optop_index_holdout(hx$models, hx$ev, hx$baseline, c = 1,
                               metrics = "deviance")
  expect_equal(ho$scores$deviance, plain$scores$deviance, tolerance = 1e-12)
})

test_that("the gain table implements Definition 1 on the actual grid", {
  hx <- optop_holdout_fixture()
  ho <- optop_index_holdout(hx$models, hx$ev, hx$baseline, c = 1)
  gt <- optop_gain_table(ho, metric = "deviance", epsilon = 0.01,
                         alpha = 0.05)

  S <- ho$scores$deviance
  z1 <- qnorm(0.95)
  for (i in 1:2) {
    delta <- S[, i + 1] - S[, i]
    delta <- delta[!is.na(delta)]
    expect_equal(gt$gains$gain[i], mean(delta), tolerance = 1e-12)
    expect_equal(gt$gains$sd[i], sd(delta), tolerance = 1e-12)
    expect_equal(gt$gains$upper[i],
                 mean(delta) + z1 * sd(delta) / sqrt(length(delta)),
                 tolerance = 1e-12)
    expect_equal(gt$gains$z[i],
                 sqrt(length(delta)) * mean(delta) / sd(delta),
                 tolerance = 1e-12)
  }
  expect_identical(gt$gains$succ_K, ho$K[-1])

  # K_hat on a fabricated elbow: gains collapse after the second grid point
  fake <- ho
  set.seed(7)
  fake$scores$deviance <- cbind(
    K2 = rnorm(10, 0.2, 0.001),
    K3 = rnorm(10, 0.5, 0.001),   # large gain 2 -> 3
    K4 = rnorm(10, 0.5, 0.001)    # negligible gain 3 -> 4
  )
  fake$K <- c(2, 3, 4)
  gt2 <- optop_gain_table(fake, metric = "deviance", epsilon = 0.01)
  expect_identical(gt2$k_hat, 3)

  expect_error(optop_gain_table(ho, metric = "nope"), "not evaluated")
  expect_error(optop_gain_table(list()), "optop_index_holdout")
})

test_that("holdout validates its inputs", {
  hx <- optop_holdout_fixture()
  bad_base <- list(pi_glob = unname(hx$baseline$pi_glob))
  expect_error(
    optop_index_holdout(hx$models, hx$ev, bad_base),
    "named"
  )
  expect_error(
    optop_index_holdout(hx$models, hx$ev, hx$baseline, conf = 1.2),
    "conf"
  )
  expect_error(
    optop_index_holdout(hx$models, hx$ev, hx$baseline, stabilize = -1),
    "stabilize"
  )
  ev_nameless <- hx$ev
  colnames(ev_nameless) <- NULL
  expect_error(
    optop_index_holdout(hx$models, ev_nameless, hx$baseline),
    "column names"
  )
})

test_that("seededlda fits fold in through the holdout layer", {
  skip_if_not_installed("seededlda")
  hx <- optop_holdout_fixture()
  tr_dfm <- quanteda::as.dfm(as.matrix(hx$tr))
  models <- lapply(2:3, function(k) {
    seededlda::textmodel_lda(tr_dfm, k = k, max_iter = 200, verbose = FALSE)
  })
  ho <- optop_index_holdout(models, hx$ev, hx$baseline, c = 1,
                            metrics = "deviance")
  expect_identical(dim(ho$scores$deviance), c(10L, 2L))
  expect_true(all(is.finite(ho$summary$macro)))
})
