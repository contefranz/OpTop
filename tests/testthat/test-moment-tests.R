# The moment-based specification tests of the paper's Section 4: instrument
# construction (training-only, exactly mean-zero, ||z||_1 = 2), the residual
# moments, and the Wald/t machinery must match hand computations.

test_that("instrument rows are mean-zero with unit-ish L1 norm", {
  hx <- optop_holdout_fixture()
  f <- as.numeric(Matrix::colSums(hx$tr))
  vocab <- colnames(hx$tr)

  z1 <- .optop_z_contrast(f, vocab)
  expect_identical(nrow(z1$Z), 1L)
  expect_equal(sum(z1$Z), 0, tolerance = 1e-14)
  expect_equal(sum(abs(z1$Z)), 2, tolerance = 1e-14)

  z2 <- .optop_z_strata(f, vocab, bins = 5)
  expect_identical(nrow(z2$Z), 4L)
  expect_equal(unname(Matrix::rowSums(z2$Z)), rep(0, 4), tolerance = 1e-14)
  expect_equal(unname(Matrix::rowSums(abs(z2$Z))), rep(2, 4),
               tolerance = 1e-14)
  # every word belongs to exactly one stratum
  expect_false(anyNA(z2$assignment))
  expect_identical(sort(unique(z2$assignment)), as.character(1:5))

  z3 <- .optop_z_fit(hx$models[[2]], hx$tr, hx$baseline, vocab,
                     strata = 4, min_doc_freq = 3, n_threads = 1L)
  expect_identical(nrow(z3$Z), 3L)
  expect_equal(unname(Matrix::rowSums(z3$Z)), rep(0, 3), tolerance = 1e-14)
  expect_equal(unname(Matrix::rowSums(abs(z3$Z))), rep(2, 3),
               tolerance = 1e-14)
  # excluded words carry zero entries in every row
  excl <- is.na(z3$assignment)
  if (any(excl)) {
    expect_true(all(abs(z3$Z[, excl]) == 0))
  }
})

test_that("the residual moments match a dense hand computation", {
  hx <- optop_holdout_fixture()
  mt <- optop_moment_test(hx$models, hx$ev, hx$tr, type = "strata",
                          bins = 4)

  # dense reconstruction of g_j = Z (d_j - theta_j phi) for one model
  i <- 2L
  m <- hx$models[[i]]
  N_al <- as.matrix(hx$ev)[, colnames(hx$tr)]
  L <- rowSums(N_al)
  stm <- slam::as.simple_triplet_matrix(N_al)
  theta <- topicmodels::posterior(m, newdata = stm)$topics
  phi <- ref_theta_phi(m)$phi
  Z <- as.matrix(mt$instruments$Z)
  G_ref <- t(apply(cbind(seq_len(nrow(N_al))), 1, function(j) {
    eps_j <- N_al[j, ] / L[j] - drop(theta[j, ] %*% phi)
    drop(Z %*% eps_j)
  }))

  expect_equal(unname(mt$moments[[i]]$g_bar), unname(colMeans(G_ref)),
               tolerance = 1e-10)
  expect_equal(unname(mt$moments[[i]]$sigma), unname(cov(G_ref)),
               tolerance = 1e-10)

  n <- nrow(G_ref)
  wald_ref <- n * drop(t(colMeans(G_ref)) %*% solve(cov(G_ref)) %*%
                         colMeans(G_ref))
  expect_equal(mt$summary$wald[i], wald_ref, tolerance = 1e-10)
  expect_equal(mt$summary$pval[i],
               pchisq(wald_ref, df = 3, lower.tail = FALSE),
               tolerance = 1e-12)

  # marginal t tests
  mar <- mt$moments[[i]]$marginal
  expect_equal(unname(mar$se), unname(sqrt(diag(cov(G_ref)) / n)),
               tolerance = 1e-12)
  expect_equal(unname(mar$t), unname(colMeans(G_ref) / mar$se),
               tolerance = 1e-10)
})

test_that("the scalar contrast satisfies W = t^2 and residuals sum to zero", {
  hx <- optop_holdout_fixture()
  mt <- optop_moment_test(hx$models[1], hx$ev, hx$tr, type = "contrast")
  expect_identical(mt$summary$q, 1L)
  expect_equal(mt$summary$wald, mt$summary$t^2, tolerance = 1e-10)

  # with d_j and p_j both probability vectors on the training support, the
  # residual sums to zero, so adding a constant to the instrument row
  # leaves the moments unchanged
  N_al <- as.matrix(hx$ev)[, colnames(hx$tr)]
  L <- rowSums(N_al)
  stm <- slam::as.simple_triplet_matrix(N_al)
  theta <- topicmodels::posterior(hx$models[[1]], newdata = stm)$topics
  phi <- ref_theta_phi(hx$models[[1]])$phi
  eps_1 <- N_al[1, ] / L[1] - drop(theta[1, ] %*% phi)
  expect_equal(sum(eps_1), 0, tolerance = 1e-10)
})

test_that("multiple-testing adjustments are applied to the marginal tests", {
  hx <- optop_holdout_fixture()
  raw <- optop_moment_test(hx$models[1], hx$ev, hx$tr, type = "strata",
                           bins = 4, adjust = "none")
  bon <- optop_moment_test(hx$models[1], hx$ev, hx$tr, type = "strata",
                           bins = 4, adjust = "bonferroni")
  m_raw <- raw$moments[[1]]$marginal
  m_bon <- bon$moments[[1]]$marginal
  expect_equal(m_raw$pval_adj, m_raw$pval, tolerance = 1e-15)
  expect_equal(m_bon$pval_adj, pmin(1, m_bon$pval * 3), tolerance = 1e-12)
})

test_that("fit-stratified moments are invariant to model-grid order", {
  hx <- optop_holdout_fixture()
  args <- list(dtm_eval = hx$ev, dtm_train = hx$tr, type = "fit",
               strata = 3, min_doc_freq = 1)
  ordered <- do.call(optop_moment_test, c(list(models = hx$models), args))
  reversed <- do.call(optop_moment_test,
                      c(list(models = rev(hx$models)), args))

  expect_identical(reversed$K, ordered$K)
  expect_equal(reversed$summary, ordered$summary, tolerance = 1e-12)
  expect_identical(names(reversed$moments), names(ordered$moments))
  expect_identical(names(reversed$instruments), names(ordered$instruments))

  for (k in as.character(ordered$K)) {
    expect_equal(reversed$moments[[k]]$g_bar, ordered$moments[[k]]$g_bar,
                 tolerance = 1e-12)
    expect_equal(reversed$moments[[k]]$sigma, ordered$moments[[k]]$sigma,
                 tolerance = 1e-12)
    expect_equal(reversed$moments[[k]]$marginal,
                 ordered$moments[[k]]$marginal, tolerance = 1e-12)
    expect_equal(as.matrix(reversed$instruments[[k]]$Z),
                 as.matrix(ordered$instruments[[k]]$Z), tolerance = 1e-12)
    expect_identical(reversed$instruments[[k]]$assignment,
                     ordered$instruments[[k]]$assignment)
  }
})

test_that("planted contamination is detected with a sane direction", {
  hx <- optop_holdout_fixture()
  # inflate the observed mass of the lowest-frequency stratum on the
  # evaluation side: the model cannot explain it, so Test 2 must reject
  f <- as.numeric(Matrix::colSums(hx$tr))
  lo_words <- order(f, seq_along(f))[1:12]
  ev <- as.matrix(hx$ev)
  ev[, lo_words] <- ev[, lo_words] + 15L
  mt <- optop_moment_test(hx$models, methods::as(ev, "CsparseMatrix"),
                          hx$tr, type = "strata", bins = 4)
  expect_true(all(mt$summary$pval < 0.01))
  # the contaminated (lowest-frequency) stratum is under-predicted:
  # positive residual moment against the high-frequency reference
  expect_true(all(vapply(mt$moments,
                         function(m) m$g_bar[1] > 0, logical(1))))
})

test_that("a singular moment covariance fails with the assumption message", {
  hx <- optop_holdout_fixture()
  # two evaluation documents cannot identify a rank-3 covariance
  expect_error(
    optop_moment_test(hx$models[1], hx$ev[1:2, ], hx$tr,
                      type = "strata", bins = 4),
    "singular"
  )
})

test_that("moment tests validate their inputs", {
  hx <- optop_holdout_fixture()
  tr_shuffled <- hx$tr[, rev(seq_len(ncol(hx$tr)))]
  expect_error(
    optop_moment_test(hx$models, hx$ev, tr_shuffled, type = "contrast"),
    "vocabulary"
  )
  expect_error(
    optop_moment_test(hx$models, hx$ev, hx$tr, type = "strata", bins = 1),
    "bins"
  )
  expect_error(
    optop_moment_test(hx$models, hx$ev, hx$tr, type = "fit", strata = 1),
    "strata"
  )
})
