# Characterization tests: the vectorized/blocked index implementations must
# agree with the naive reference implementation (helper-reference.R) that
# follows the paper's definitions with explicit loops.

optop_index_fun <- function(metric) {
  switch(metric,
         se = optop_index_se,
         chisq = optop_index_chisq,
         deviance = optop_index_deviance)
}

test_that("document-level indices match the naive reference", {
  withr::local_options(lifecycle_verbosity = "quiet")
  fx <- optop_test_fixture()
  for (metric in c("se", "chisq", "deviance")) {
    fun <- optop_index_fun(metric)
    for (m in fx$models) {
      got <- fun(m, fx$dtm, fx$partition, fx$baseline,
                 macro = TRUE, ztest = TRUE)
      ref <- ref_index_document(m, fx$dtm, fx$partition, fx$baseline, metric)

      info <- sprintf("metric=%s K=%d", metric, ref$K)
      expect_equal(got$r2, ref$r2, tolerance = 1e-10, info = info)
      expect_equal(got$r2_macro, ref$r2_macro, tolerance = 1e-10, info = info)
      expect_equal(got$r2_doc, ref$r2_doc, tolerance = 1e-10, info = info)
      expect_identical(got$K, ref$K)
      expect_identical(got$metric, metric)

      ref_z <- ref_ztest(ref$r2_doc[ref$D_null > 0])
      expect_equal(got$ztest$z, ref_z$z, tolerance = 1e-10, info = info)
      expect_equal(got$ztest$pval, ref_z$pval, tolerance = 1e-10, info = info)
      expect_equal(got$ztest$se, ref_z$se, tolerance = 1e-10, info = info)
      expect_equal(got$ztest$ci, ref_z$ci, tolerance = 1e-10, info = info)
      expect_identical(got$ztest$J, ref_z$J)
    }
  }
})

test_that("word-level indices match the naive reference", {
  fx <- optop_test_fixture()
  for (metric in c("se", "chisq", "deviance")) {
    fun <- optop_index_fun(metric)
    for (m in fx$models) {
      got <- fun(m, fx$dtm, fx$partition, fx$baseline, level = "word")
      ref <- ref_index_word(m, fx$dtm, fx$partition, fx$baseline, metric)

      info <- sprintf("metric=%s K=%d", metric, ref$K)
      expect_equal(got$r2_word, ref$r2_word, tolerance = 1e-10, info = info)
      expect_equal(got$r2_micro_word, ref$r2_micro_word,
                   tolerance = 1e-10, info = info)
      expect_equal(got$r2_macro_word, ref$r2_macro_word,
                   tolerance = 1e-10, info = info)
      expect_identical(got$K, ref$K)
    }
  }
})

test_that("results are invariant to block_size at both levels", {
  fx <- optop_test_fixture()
  m <- fx$models[[2]]
  for (metric in c("se", "chisq", "deviance")) {
    fun <- optop_index_fun(metric)

    doc_default <- fun(m, fx$dtm, fx$partition, fx$baseline, macro = TRUE)
    word_default <- fun(m, fx$dtm, fx$partition, fx$baseline, level = "word")
    word_small <- fun(m, fx$dtm, fx$partition, fx$baseline,
                      level = "word", block_size = 13L)

    expect_equal(word_small$r2_word, word_default$r2_word, tolerance = 1e-12)
    expect_equal(word_small$r2_micro_word, word_default$r2_micro_word,
                 tolerance = 1e-12)
    expect_equal(word_small$r2_macro_word, word_default$r2_macro_word,
                 tolerance = 1e-12)

    # document-level block size is derived internally; re-run to confirm
    # determinism of the blocked path
    doc_again <- fun(m, fx$dtm, fx$partition, fx$baseline, macro = TRUE)
    expect_identical(doc_again$r2, doc_default$r2)
    expect_identical(doc_again$r2_doc, doc_default$r2_doc)
  }
})

test_that("baseline-topic augmentation is inert when reopt = 'none'", {
  withr::local_options(lifecycle_verbosity = "quiet")
  fx <- optop_test_fixture()
  m <- fx$models[[1]]
  for (metric in c("se", "chisq")) {
    fun <- optop_index_fun(metric)
    with_aug <- fun(m, fx$dtm, fx$partition, fx$baseline,
                    macro = TRUE, add_baseline_topic = TRUE)
    without_aug <- fun(m, fx$dtm, fx$partition, fx$baseline,
                       macro = TRUE, add_baseline_topic = FALSE)
    expect_equal(with_aug$r2, without_aug$r2, tolerance = 1e-12)
    expect_equal(with_aug$r2_doc, without_aug$r2_doc, tolerance = 1e-12)
    expect_equal(with_aug$r2_macro, without_aug$r2_macro, tolerance = 1e-12)
  }
})

test_that("SE re-optimization matches the reference and bounds R2 below by 0", {
  withr::local_options(lifecycle_verbosity = "quiet")
  fx <- optop_test_fixture()
  for (m in fx$models) {
    got <- optop_index_se(m, fx$dtm, fx$partition, fx$baseline,
                          macro = TRUE, reopt = "se")
    ref <- ref_index_document(m, fx$dtm, fx$partition, fx$baseline,
                              metric = "se", reopt = "se")
    expect_equal(got$r2, ref$r2, tolerance = 1e-10)
    expect_equal(got$r2_doc, ref$r2_doc, tolerance = 1e-10)
    expect_equal(got$r2_macro, ref$r2_macro, tolerance = 1e-10)
    # lambda in [0,1] guarantees per-document non-negativity
    expect_true(all(got$r2_doc >= -1e-12))
    expect_gte(got$r2, -1e-12)
  }
})

test_that("macro and ztest components are only returned when requested", {
  withr::local_options(lifecycle_verbosity = "quiet")
  fx <- optop_test_fixture()
  m <- fx$models[[1]]
  res <- optop_index_se(m, fx$dtm, fx$partition, fx$baseline)
  expect_null(res$r2_macro)
  expect_null(res$ztest)

  res_macro <- optop_index_se(m, fx$dtm, fx$partition, fx$baseline,
                              macro = TRUE, ztest = TRUE)
  expect_true(is.numeric(res_macro$r2_macro))
  expect_named(res_macro$ztest, c("z", "pval", "se", "ci", "J"))
  expect_true(res_macro$ztest$pval >= 0 && res_macro$ztest$pval <= 1)
})

test_that("precomputed null discrepancies reproduce the self-computed ones", {
  fx <- optop_test_fixture()
  m <- fx$models[[2]]
  pi_row <- fx$baseline$pi_glob[colnames(fx$dtm)]
  impls <- list(se = OpTop:::.optop_index_se_impl,
                chisq = OpTop:::.optop_index_chisq_impl,
                deviance = OpTop:::.optop_index_deviance_impl)

  for (metric in names(impls)) {
    impl <- impls[[metric]]
    args_common <- if (metric == "deviance") {
      list(m, fx$dtm, fx$partition, fx$baseline, macro = TRUE, reopt = "none")
    } else {
      list(m, fx$dtm, fx$partition, fx$baseline, macro = TRUE, reopt = "none",
           add_baseline_topic = TRUE)
    }

    for (lvl in c("document", "word")) {
      null_disc <- OpTop:::.optop_index_null(fx$dtm, fx$partition, pi_row,
                                             metric = metric, level = lvl)
      plain <- do.call(impl, c(args_common, list(level = lvl, block_size = NULL,
                                                 ztest = FALSE)))
      hoisted <- do.call(impl, c(args_common, list(level = lvl, block_size = NULL,
                                                   ztest = FALSE,
                                                   null_disc = null_disc)))
      expect_identical(hoisted, plain,
                       info = sprintf("metric=%s level=%s", metric, lvl))
    }
  }
})

test_that("misaligned inputs raise the documented errors", {
  fx <- optop_test_fixture()
  m <- fx$models[[1]]

  # shuffled dtm columns
  dtm_shuffled <- fx$dtm[, rev(seq_len(fx$W))]
  expect_error(
    optop_index_se(m, dtm_shuffled, fx$partition, fx$baseline),
    "DTM vocabulary/order differs from model"
  )

  # baseline with wrong vocabulary
  bad_baseline <- fx$baseline
  names(bad_baseline$pi_glob) <- paste0("x_", names(bad_baseline$pi_glob))
  expect_error(
    optop_index_se(m, fx$dtm, fx$partition, bad_baseline),
    "Baseline vocabulary does not match model"
  )

  # unnamed baseline of the wrong length
  short_baseline <- list(pi_glob = unname(fx$baseline$pi_glob[-1]))
  expect_error(
    optop_index_se(m, fx$dtm, fx$partition, short_baseline),
    "Baseline length != model vocab"
  )

  # partition missing vocabulary labels
  bad_partition <- fx$partition
  colnames(bad_partition$rare_mask) <- NULL
  expect_error(
    optop_index_se(m, fx$dtm, bad_partition, fx$baseline),
    "Partition vocabulary != model vocabulary"
  )
})

test_that("the index engine kernels validate their inputs", {
  J <- 3L; W <- 4L
  N <- Matrix::Matrix(matrix(1, J, W), sparse = TRUE)
  N_t <- Matrix::t(N)
  E <- matrix(1 / W, J, W)
  L <- rep(4, J)
  pi_w <- rep(1 / W, W)

  # word-level kernel
  expect_error(
    optop_index_word_core(E[-1, ], N, 0L, W, L, pi_w, 1e-8,
                          TRUE, TRUE, TRUE, FALSE, FALSE, 1L),
    "E_block"
  )
  expect_error(
    optop_index_word_core(E, N, 0L, W, L[-1], pi_w, 1e-8,
                          TRUE, TRUE, TRUE, FALSE, FALSE, 1L),
    "L must have one entry per document"
  )
  expect_error(
    optop_index_word_core(E, N, 0L, W, L, pi_w[-1], 1e-8,
                          TRUE, TRUE, TRUE, FALSE, FALSE, 1L),
    "pi_w"
  )

  # document-level kernel
  tww <- matrix(1 / W, W, 2L)
  theta_blk <- matrix(1 / 2, J, 2L)
  bits <- as.raw(rep(0L, J * W))
  ok <- rep(TRUE, J)
  expect_error(
    optop_index_doc_core(tww[-1, ], theta_blk, N_t, 0L, bits, L, pi_w, ok, 1e-8,
                         TRUE, TRUE, TRUE, FALSE, FALSE, 1L),
    "one row per feature"
  )
  expect_error(
    optop_index_doc_core(tww, theta_blk[-1, ], N_t, 0L, bits, L, pi_w, ok, 1e-8,
                         TRUE, TRUE, TRUE, FALSE, FALSE, 1L),
    "theta_blk"
  )
  expect_error(
    optop_index_doc_core(tww, theta_blk, N_t, 0L, bits, L, pi_w[-1], ok, 1e-8,
                         TRUE, TRUE, TRUE, FALSE, FALSE, 1L),
    "pi_row"
  )
  expect_error(
    optop_index_doc_core(tww, theta_blk, N_t, 0L, bits, L, pi_w, ok[-1], 1e-8,
                         TRUE, TRUE, TRUE, FALSE, FALSE, 1L),
    "chisq_min_ok"
  )
})
