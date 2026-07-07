test_that("optop_make_partition returns lengths, mask, and the Pearson rule", {
  fx <- optop_test_fixture()
  part <- fx$partition

  expect_named(part, c("rare_mask", "L", "chisq_min_ok",
                       "chisq_min_report", "c"))
  expect_identical(dim(part$rare_mask), c(fx$J, fx$W))
  expect_type(part$rare_mask, "logical")
  expect_identical(colnames(part$rare_mask), colnames(fx$counts))
  expect_identical(rownames(part$rare_mask), rownames(fx$counts))
  expect_equal(unname(part$L), unname(Matrix::rowSums(fx$dtm)))
  expect_type(part$chisq_min_ok, "logical")
  expect_length(part$chisq_min_ok, fx$J)
  expect_named(part$chisq_min_report,
               c("n_excluded", "share", "excluded_mass"))
  expect_identical(part$c, 5)

  # the fixture must exercise both branches of the harmonized support
  expect_true(any(part$rare_mask))
  expect_true(any(!part$rare_mask))
})

test_that("rare_mask is the harmonized union of baseline and model grid", {
  fx <- optop_test_fixture()

  L <- Matrix::rowSums(fx$dtm)
  tau <- 5 / pmax(L, 1)
  pi_glob <- as.numeric(Matrix::colSums(fx$dtm)) / sum(fx$dtm)
  i_min <- Reduce(pmin, lapply(fx$models, function(m) {
    tp <- ref_theta_phi(m)
    tp$theta %*% tp$phi
  }))
  # the no-topics baseline enters the union: rare iff
  # min(pi_glob(w), min_K i^K_jw) < tau_j
  i_min <- pmin(i_min, matrix(pi_glob, nrow = fx$J, ncol = fx$W, byrow = TRUE))
  expected <- i_min < tau  # recycles tau_j down each column (row j)

  expect_equal(unname(fx$partition$rare_mask), unname(expected))
})

test_that("chisq_min_ok equals the naive grid-wide inclusion rule", {
  fx <- optop_test_fixture()
  part <- fx$partition
  L <- part$L
  pi_glob <- as.numeric(Matrix::colSums(fx$dtm)) / sum(fx$dtm)

  E_min_by_model <- sapply(fx$models, function(m) {
    tp <- ref_theta_phi(m)
    I <- tp$theta %*% tp$phi
    L * rowSums(I * part$rare_mask)
  })
  E_min_min <- apply(E_min_by_model, 1, min)
  B_min <- L * as.numeric(part$rare_mask %*% pi_glob)

  expect_equal(unname(part$chisq_min_ok),
               unname(pmin(E_min_min, B_min) >= 5))
  has_min <- rowSums(part$rare_mask) > 0
  excluded <- has_min & !part$chisq_min_ok
  expect_identical(part$chisq_min_report$n_excluded, sum(excluded))
  expect_equal(part$chisq_min_report$share, mean(excluded))
})

test_that("partition is invariant to the block size", {
  fx <- optop_test_fixture()
  part_small <- optop_make_partition(fx$models, fx$dtm, c = 5, block = 7)
  expect_equal(part_small$rare_mask, fx$partition$rare_mask)
  expect_equal(part_small$L, fx$partition$L)
  expect_identical(part_small$chisq_min_ok, fx$partition$chisq_min_ok)
})

test_that("larger c marks more words as rare", {
  fx <- optop_test_fixture()
  part_loose <- optop_make_partition(fx$models, fx$dtm, c = 20)
  expect_true(all(fx$partition$rare_mask <= part_loose$rare_mask))
  expect_gt(sum(part_loose$rare_mask), sum(fx$partition$rare_mask))
})

test_that("the partition kernels validate their inputs and clamp the threads", {
  J <- 3L; W <- 4L
  mask <- matrix(FALSE, J, W)
  theta <- matrix(1 / 2, J, 2L)
  phi <- matrix(1 / W, 2L, W)
  pi_glob <- rep(1 / W, W)
  tau <- rep(0.1, J)

  expect_error(
    optop_partition_fill_core(mask, list(theta), list(), pi_glob, tau, 2L, 1L),
    "one entry per model"
  )
  expect_error(
    optop_partition_fill_core(mask, list(theta), list(phi), pi_glob[-1],
                              tau, 2L, 1L),
    "pi_glob"
  )
  expect_error(
    optop_partition_fill_core(mask, list(theta), list(phi), pi_glob,
                              tau[-1], 2L, 1L),
    "tau"
  )
  expect_error(
    optop_partition_fill_core(mask, list(theta), list(phi), pi_glob,
                              tau, 0L, 1L),
    "block"
  )
  expect_error(
    optop_partition_fill_core(mask, list(theta[-1, ]), list(phi), pi_glob,
                              tau, 2L, 1L),
    "one row per document"
  )
  expect_error(
    optop_partition_fill_core(mask, list(theta), list(phi[, -1]), pi_glob,
                              tau, 2L, 1L),
    "one column per feature"
  )
  expect_error(
    optop_partition_fill_core(mask, list(theta), list(rbind(phi, phi[1, ])),
                              pi_glob, tau, 2L, 1L),
    "number of topics"
  )

  # a non-positive thread count is clamped to one thread (in-place fill)
  mask0 <- matrix(FALSE, J, W)
  mask1 <- matrix(FALSE, J, W)
  optop_partition_fill_core(mask0, list(theta), list(phi), pi_glob, tau, 2L, 0L)
  optop_partition_fill_core(mask1, list(theta), list(phi), pi_glob, tau, 2L, 1L)
  expect_identical(mask0, mask1)

  # rare-mass kernel guards
  lmask <- matrix(TRUE, J, W)
  expect_error(
    optop_partition_minmass_core(lmask, list(theta), list(), 2L, 1L),
    "one entry per model"
  )
  expect_error(
    optop_partition_minmass_core(lmask, list(theta), list(phi), 0L, 1L),
    "block"
  )
  expect_error(
    optop_partition_minmass_core(lmask, list(theta[-1, ]), list(phi), 2L, 1L),
    "one row per document"
  )
  expect_error(
    optop_partition_minmass_core(lmask, list(theta), list(phi[, -1]), 2L, 1L),
    "one column per feature"
  )
  expect_error(
    optop_partition_minmass_core(lmask, list(theta),
                                 list(rbind(phi, phi[1, ])), 2L, 1L),
    "number of topics"
  )

  # with an all-rare mask the rare mass is the full row sum of theta phi
  got <- optop_partition_minmass_core(lmask, list(theta), list(phi), 2L, 0L)
  expect_equal(unname(got[, 1]), unname(rowSums(theta %*% phi)))
})
