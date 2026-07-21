test_that("optop_make_partition returns the sparse format and the Pearson rule", {
  fx <- optop_test_fixture()
  part <- fx$partition

  expect_named(part, c("nonrare_offsets", "nonrare_words", "vocab", "L",
                       "chisq_min_ok", "chisq_min_report", "c", "format"))
  expect_identical(part$format, 2L)
  expect_type(part$nonrare_offsets, "double")
  expect_length(part$nonrare_offsets, fx$J + 1L)
  expect_identical(part$nonrare_offsets[1], 0)
  expect_true(all(diff(part$nonrare_offsets) >= 0))
  expect_type(part$nonrare_words, "integer")
  expect_length(part$nonrare_words, as.integer(part$nonrare_offsets[fx$J + 1]))
  expect_identical(part$vocab, colnames(fx$counts))
  expect_equal(unname(part$L), unname(Matrix::rowSums(fx$dtm)))
  expect_type(part$chisq_min_ok, "logical")
  expect_length(part$chisq_min_ok, fx$J)
  expect_named(part$chisq_min_report,
               c("n_excluded", "share", "excluded_mass"))
  expect_identical(part$c, 5)

  # word indices are 0-based, in range, ascending within each document
  expect_true(all(part$nonrare_words >= 0L))
  expect_true(all(part$nonrare_words < fx$W))
  for (j in seq_len(fx$J)) {
    o0 <- part$nonrare_offsets[j]
    o1 <- part$nonrare_offsets[j + 1]
    if (o1 - o0 > 1) {
      expect_false(is.unsorted(part$nonrare_words[(o0 + 1):o1],
                               strictly = TRUE))
    }
  }

  # the fixture must exercise both branches of the harmonized support
  mask <- ref_rare_mask(part, fx$J, fx$W)
  expect_true(any(mask))
  expect_true(any(!mask))
})

test_that("the rare set is the harmonized union of baseline and model grid", {
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

  got <- ref_rare_mask(fx$partition, fx$J, fx$W)
  expect_equal(unname(got), unname(expected))
})

test_that("chisq_min_ok equals the naive grid-wide inclusion rule", {
  fx <- optop_test_fixture()
  part <- fx$partition
  L <- part$L
  pi_glob <- as.numeric(Matrix::colSums(fx$dtm)) / sum(fx$dtm)
  mask <- ref_rare_mask(part, fx$J, fx$W)

  E_min_by_model <- sapply(fx$models, function(m) {
    tp <- ref_theta_phi(m)
    I <- tp$theta %*% tp$phi
    L * rowSums(I * mask)
  })
  E_min_min <- apply(E_min_by_model, 1, min)
  B_min <- L * as.numeric(mask %*% pi_glob)

  expect_identical(unname(part$chisq_min_ok),
                   unname(pmin(E_min_min, B_min) >= 5))
  has_min <- rowSums(mask) > 0
  excluded <- has_min & !part$chisq_min_ok
  expect_identical(part$chisq_min_report$n_excluded, sum(excluded))
  expect_equal(part$chisq_min_report$share, mean(excluded))
})

test_that("the deprecated block argument is accepted and inert", {
  fx <- optop_test_fixture()
  part_small <- optop_make_partition(fx$models, fx$dtm, c = 5, block = 7)
  expect_equal(part_small$nonrare_offsets, fx$partition$nonrare_offsets)
  expect_identical(part_small$nonrare_words, fx$partition$nonrare_words)
  expect_equal(part_small$L, fx$partition$L)
  expect_identical(part_small$chisq_min_ok, fx$partition$chisq_min_ok)
})

test_that("larger c marks more words as rare", {
  fx <- optop_test_fixture()
  part_loose <- optop_make_partition(fx$models, fx$dtm, c = 20)
  mask_tight <- ref_rare_mask(fx$partition, fx$J, fx$W)
  mask_loose <- ref_rare_mask(part_loose, fx$J, fx$W)
  expect_true(all(mask_tight <= mask_loose))
  expect_gt(sum(mask_loose), sum(mask_tight))
})

test_that("candidate prefixes match brute force over random corpora", {
  # the exactness of the prefilter: non-rare requires pi_glob(w) >= tau_j
  # (the null belongs to the augmented union), so the sparse construction
  # must reproduce the brute-force dense rule cell by cell
  set.seed(202)
  for (rep in 1:5) {
    J <- sample(5:20, 1)
    W <- sample(10:60, 1)
    K <- sample(2:4, 1)
    counts <- matrix(rpois(J * W, 2), J, W,
                     dimnames = list(sprintf("d%02d", 1:J),
                                     sprintf("w%02d", 1:W)))
    counts[1, ] <- counts[1, ] + 1  # no empty docs
    dtm <- methods::as(Matrix::Matrix(counts, sparse = TRUE),
                       "CsparseMatrix")
    rdir <- function(n, k) {
      g <- matrix(stats::rgamma(n * k, 1), n)
      g / rowSums(g)
    }
    models <- lapply(seq_len(2), function(m) {
      theta <- rdir(J, K + m)
      dimnames(theta) <- list(rownames(counts), NULL)
      phi <- rdir(K + m, W)
      colnames(phi) <- colnames(counts)
      OpTop:::.optop_tp(theta, phi)
    })
    cc <- sample(c(1, 2, 5), 1)
    part <- optop_make_partition(models, dtm, c = cc)

    L <- Matrix::rowSums(dtm)
    tau <- cc / pmax(L, 1)
    pi_glob <- as.numeric(Matrix::colSums(dtm)) / sum(dtm)
    i_min <- Reduce(pmin, lapply(models, function(m) m$theta %*% m$phi))
    i_min <- pmin(i_min, matrix(pi_glob, J, W, byrow = TRUE))
    expected <- i_min < tau

    expect_equal(unname(ref_rare_mask(part, J, W)), unname(expected))
  }
})

test_that("dense-mask partitions upgrade automatically with an alert", {
  fx <- optop_test_fixture()
  mask <- ref_rare_mask(fx$partition, fx$J, fx$W)
  legacy <- list(rare_mask = mask, L = fx$partition$L,
                 chisq_min_ok = fx$partition$chisq_min_ok,
                 chisq_min_report = fx$partition$chisq_min_report,
                 c = fx$partition$c)

  expect_message(
    up <- OpTop:::.optop_partition_upgrade(legacy),
    "sparse format"
  )
  expect_identical(up$format, 2L)
  expect_equal(up$nonrare_offsets, fx$partition$nonrare_offsets)
  expect_identical(up$nonrare_words, fx$partition$nonrare_words)
  expect_null(up$rare_mask)
  expect_identical(up$vocab, fx$partition$vocab)

  # the upgraded partition reproduces the index of the native one exactly
  base <- optop_make_baseline(fx$dtm)
  res_native <- suppressMessages(
    optop_index_deviance(fx$models[[2]], fx$dtm, fx$partition, base,
                         macro = TRUE)
  )
  res_legacy <- suppressMessages(
    optop_index_deviance(fx$models[[2]], fx$dtm, legacy, base, macro = TRUE)
  )
  expect_identical(res_native$r2, res_legacy$r2)
  expect_identical(res_native$r2_doc, res_legacy$r2_doc)
})

test_that("the partition kernels validate their inputs", {
  J <- 3L; W <- 4L
  theta <- matrix(1 / 2, J, 2L)
  phi <- matrix(1 / W, 2L, W)
  tau <- rep(0.1, J)
  ord <- as.integer(seq_len(W) - 1L)
  cand_off <- as.numeric(c(0, W, 2 * W, 3 * W))
  keep <- as.raw(rep(1L, 3 * W))

  expect_error(
    optop_partition_pass_core(theta, rbind(phi, phi[1, ]), ord, cand_off,
                              keep, tau, 1L),
    "number of topics"
  )
  expect_error(
    optop_partition_pass_core(theta[-1, ], phi, ord, cand_off, keep, tau, 1L),
    "one row per document"
  )
  expect_error(
    optop_partition_pass_core(theta, phi, ord, cand_off[-1], keep, tau, 1L),
    "J \\+ 1"
  )
  expect_error(
    optop_partition_sums_core(theta[-1, ], phi,
                              as.numeric(c(0, 1, 2, 3)),
                              as.integer(c(0, 1, 2)), 1L),
    "one row per document"
  )

  # thread clamping: a non-positive count equals one thread
  a <- optop_partition_candidates_core(sort(rep(1 / W, W), decreasing = TRUE),
                                       tau, 0L)
  b <- optop_partition_candidates_core(sort(rep(1 / W, W), decreasing = TRUE),
                                       tau, 1L)
  expect_identical(a, b)
})
