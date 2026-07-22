# Characterization tests for optop_select(): the C++ core must agree with the
# naive per-document reference (helper-reference.R) for the Test 1 statistic
# of Eq. (8) consumed by the "sequential" and "min" rules.

# numeric tests run silent (the wprop fixture lives in helper-fixtures.R):
# the unconditional document-drop alert is muffled where it is not the
# behavior under test
run_optimal_topic <- function(...) {
  suppressMessages(optop_select(..., do_plot = FALSE))
}

expect_matches_reference <- function(got, ref, info = NULL) {
  expect_s3_class(got, "data.table")
  expect_identical(names(got), c("topic", "OpTop", "df", "pval"))
  expect_equal(got$topic, ref$topic, info = info)
  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10, info = info)
  expect_equal(got$df, ref$df, info = info)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10, info = info)
}

test_that("optop_select() matches the naive Eq. (8) reference", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  for (q in c(0.95, 0.80)) {
    got <- run_optimal_topic(fx$models, wp$wdfm, q = q)
    ref <- ref_optimal_topic(fx$models, wp$W_prop, q = q)
    expect_matches_reference(got, ref, info = sprintf("q=%.2f", q))
  }
})

test_that("the returned table does not depend on the selection rule", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  got_seq <- run_optimal_topic(fx$models, wp$wdfm)
  got_min <- run_optimal_topic(fx$models, wp$wdfm, selection = "min")
  expect_identical(got_seq, got_min)
})

test_that("selection = \"legacy\" was removed in 0.16.0", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_error(
    run_optimal_topic(fx$models, wp$wdfm, selection = "legacy"),
    "should be one of"
  )
})

test_that("the selection rules pick the expected K", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  ref <- ref_optimal_topic(fx$models, wp$W_prop, q = 0.95)

  # sequential: smallest K the test fails to reject
  seq_pick <- ref$topic[ref$pval > 0.05][1]
  expect_false(is.na(seq_pick))
  expect_message(
    optop_select(fx$models, wp$wdfm, do_plot = FALSE, verbose = TRUE),
    sprintf("Optimal model has %d topics", seq_pick)
  )

  # min: global minimum of the standardized statistic
  min_pick <- ref$topic[which.min(ref$OpTop)]
  expect_message(
    optop_select(fx$models, wp$wdfm, selection = "min",
                 do_plot = FALSE, verbose = TRUE),
    sprintf("Optimal model has %d topics", min_pick)
  )

  # alpha = 1 makes 'pval > alpha' unsatisfiable: the sequential rule must
  # warn and fall back to the global minimum
  expect_message(
    optop_select(fx$models, wp$wdfm, alpha = 1, do_plot = FALSE),
    "falling back to the global minimum"
  )
  expect_message(
    optop_select(fx$models, wp$wdfm, alpha = 1, do_plot = FALSE,
                 verbose = TRUE),
    sprintf("Optimal model has %d topics", min_pick)
  )
})

test_that("documents missing from the models are dropped, order preserved", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  counts_extra <- rbind(fx$counts, docXX = fx$counts[1, ])
  wdfm_extra <- quanteda::dfm_weight(quanteda::as.dfm(counts_extra),
                                     scheme = "prop")

  got <- run_optimal_topic(fx$models, wdfm_extra)
  ref <- ref_optimal_topic(fx$models, wp$W_prop, q = 0.95)

  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10)
})

test_that("optop_select() is invariant to document order in the dfm", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  set.seed(42)
  perm <- sample(nrow(fx$counts))
  wdfm_perm <- quanteda::dfm_weight(quanteda::as.dfm(fx$counts[perm, ]),
                                    scheme = "prop")

  got <- run_optimal_topic(fx$models, wdfm_perm)
  ref <- ref_optimal_topic(fx$models, wp$W_prop, q = 0.95)

  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10)
})

test_that("optop_select() is chatty by default and silent with verbose = FALSE", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_no_message(optop_select(fx$models, wp$wdfm, do_plot = FALSE,
                                  verbose = FALSE))
  expect_message(
    optop_select(fx$models, wp$wdfm, do_plot = FALSE),
    "Optimal model has"
  )
})

test_that("dropping unmatched documents is signalled even when silent", {
  fx <- optop_test_fixture()

  counts_extra <- rbind(fx$counts, docXX = fx$counts[1, ])
  wdfm_extra <- quanteda::dfm_weight(quanteda::as.dfm(counts_extra),
                                     scheme = "prop")

  expect_message(
    optop_select(fx$models, wdfm_extra, do_plot = FALSE),
    "Removed 1 document"
  )
})

test_that("optimal_topic() delegates to optop_select() with a one-time note", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  ref <- optop_select(fx$models, wp$wdfm, do_plot = FALSE, verbose = FALSE)

  # clear the session flag so the note is observable in this test
  rm(list = ls(OpTop:::.optop_rename_env), envir = OpTop:::.optop_rename_env)
  expect_message(
    got <- optimal_topic(fx$models, wp$wdfm, do_plot = FALSE,
                         verbose = FALSE),
    "renamed optop_select"
  )
  expect_identical(got, ref)

  # the note fires once per session; later calls are silent
  expect_silent(
    optimal_topic(fx$models, wp$wdfm, do_plot = FALSE, verbose = FALSE)
  )
})

test_that("the optop.cache_mb option bounds the cache without changing results", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  ref <- run_optimal_topic(fx$models, wp$wdfm)

  # a zero budget disables caching: every model is re-extracted, results
  # are identical because extraction is deterministic
  old <- options(optop.cache_mb = 0)
  on.exit(options(old), add = TRUE)
  got <- run_optimal_topic(fx$models, wp$wdfm)
  expect_identical(got, ref)

  options(optop.cache_mb = "big")
  expect_error(run_optimal_topic(fx$models, wp$wdfm), "optop.cache_mb")
})

test_that("optop_select() validates its inputs", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_error(optop_select("not a list", wp$wdfm),
               "must be a list")
  expect_error(optop_select(fx$models[1L], wp$wdfm),
               "multiple topic models")
  expect_error(optop_select(list(1, 2), wp$wdfm),
               "unsupported topic model class")
  expect_error(optop_select(fx$models, fx$counts),
               "must be a dfm")
  expect_error(optop_select(fx$models, wp$wdfm, q = "a"),
               "q must be a numeric")
  expect_error(optop_select(fx$models, wp$wdfm, alpha = "a"),
               "alpha must be a numeric")
  expect_error(optop_select(fx$models, wp$wdfm, selection = "typo"),
               "should be one of")
})
