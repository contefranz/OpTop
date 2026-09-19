# The provenance guard of the word-level deviance null (0.20.1): the
# compiled kernel omits the linear term of the null side, an omission
# that is exact only in sample. The guard compares the baseline against
# the pooled column sums and warns (class "optop_word_null_baseline")
# when the in-sample precondition fails; it never touches the numeric
# path, so every in-sample result stays bit-identical.

test_that("the guard predicate accepts in-sample witnesses and rejects foreign ones", {
  fx <- optop_test_fixture()
  cs <- as.numeric(Matrix::colSums(fx$dtm))

  expect_no_warning(
    ok <- OpTop:::.optop_check_word_null_baseline(fx$dtm, fx$partition$L,
                                                  fx$baseline$pi_glob)
  )
  expect_true(ok)

  # precomputed column sums and a scalar length total are equivalent
  expect_no_warning(
    ok2 <- OpTop:::.optop_check_word_null_baseline(cs, sum(fx$partition$L),
                                                   fx$baseline$pi_glob)
  )
  expect_true(ok2)

  # a corpus input pools its shards
  corp <- optop_corpus(list(fx$dtm[1:12, ], fx$dtm[13:30, ]))
  expect_no_warning(
    ok3 <- OpTop:::.optop_check_word_null_baseline(corp, fx$partition$L,
                                                   fx$baseline$pi_glob)
  )
  expect_true(ok3)

  # a perturbed baseline fails with the stable class
  pi_bad <- fx$baseline$pi_glob
  pi_bad[1] <- pi_bad[1] * 1.5
  pi_bad <- pi_bad / sum(pi_bad)
  expect_warning(
    bad <- OpTop:::.optop_check_word_null_baseline(fx$dtm, fx$partition$L,
                                                   pi_bad),
    class = "optop_word_null_baseline"
  )
  expect_false(bad)

  # a structurally foreign baseline (wrong length) fails too
  expect_warning(
    OpTop:::.optop_check_word_null_baseline(fx$dtm, fx$partition$L,
                                            pi_bad[-1]),
    class = "optop_word_null_baseline"
  )
})

test_that("in-sample word-level deviance is silent and matches the reference", {
  fx <- optop_test_fixture()
  m <- fx$models[[2]]

  expect_no_warning(
    res <- optop_index_deviance(m, fx$dtm, fx$partition, fx$baseline,
                                level = "word")
  )
  ref <- ref_index_word(m, fx$dtm, fx$partition, fx$baseline, "deviance")
  expect_equal(res$r2_word, ref$r2_word, tolerance = 1e-10)
  expect_equal(res$r2_micro_word, ref$r2_micro_word, tolerance = 1e-10)

  # sharded in-sample evaluation regroups sums; the tolerance absorbs it
  corp <- optop_corpus(list(fx$dtm[1:12, ], fx$dtm[13:30, ]))
  expect_no_warning(
    optop_index_deviance(m, corp, fx$partition, fx$baseline, level = "word")
  )
})

test_that("foreign and smoothed baselines warn at the word level for deviance", {
  fx <- optop_test_fixture()
  m <- fx$models[[2]]

  # a baseline from a subset of the corpus: same vocabulary, different
  # distribution; the call is answered, with the warning
  base_sub <- optop_make_baseline(fx$dtm[1:15, ])
  expect_warning(
    res <- optop_index_deviance(m, fx$dtm, fx$partition, base_sub,
                                level = "word"),
    class = "optop_word_null_baseline"
  )
  expect_true(all(is.finite(res$d_model)))

  # smoothing breaks the same identity
  base_smooth <- optop_make_baseline(fx$dtm, smooth_lambda = 0.5)
  expect_warning(
    optop_index_deviance(m, fx$dtm, fx$partition, base_smooth,
                         level = "word"),
    class = "optop_word_null_baseline"
  )
})

test_that("the guard is scoped to the word-level deviance null", {
  fx <- optop_test_fixture()
  m <- fx$models[[2]]
  base_smooth <- optop_make_baseline(fx$dtm, smooth_lambda = 0.5)

  # the squared-error and Pearson word-level nulls carry no linear term
  expect_no_warning(
    optop_index_se(m, fx$dtm, fx$partition, base_smooth, level = "word")
  )
  expect_no_warning(
    optop_index_chisq(m, fx$dtm, fx$partition, base_smooth, level = "word")
  )

  # the document level is immune by construction (min-bin complement)
  expect_no_warning(
    suppressMessages(
      optop_index_deviance(m, fx$dtm, fx$partition, base_smooth,
                           macro = TRUE)
    )
  )
})

test_that("the grid table warns once per call, not per model", {
  fx <- optop_test_fixture()
  base_smooth <- optop_make_baseline(fx$dtm, smooth_lambda = 0.5)

  cw <- testthat::capture_warnings(
    optop_index_table(fx$models[1:2], fx$dtm, metrics = "deviance",
                      level = "word", partition = fx$partition,
                      baseline = base_smooth)
  )
  expect_length(cw, 1L)

  expect_no_warning(
    optop_index_table(fx$models[1:2], fx$dtm, metrics = "deviance",
                      level = "word", partition = fx$partition,
                      baseline = fx$baseline)
  )
})

test_that("the moment test guards its fit-stratified score", {
  hx <- optop_holdout_fixture()
  m1 <- hx$models[1]

  expect_no_warning(
    suppressMessages(
      optop_moment_test(m1, hx$ev, hx$tr, type = "fit")
    )
  )

  base_smooth <- optop_make_baseline(hx$tr, smooth_lambda = 0.5)
  cw <- testthat::capture_warnings(
    suppressMessages(
      optop_moment_test(hx$models, hx$ev, hx$tr, type = "fit",
                        baseline = base_smooth)
    )
  )
  expect_length(cw[grepl("word-level deviance null", cw)], 1L)

  # the strata instruments never read the word-level null
  expect_no_warning(
    suppressMessages(
      optop_moment_test(m1, hx$ev, hx$tr, type = "strata",
                        baseline = base_smooth)
    )
  )
})
