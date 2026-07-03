# Characterization tests for optimal_topic(): the C++ core must agree with the
# naive per-document references (helper-reference.R) — the calibrated Test 1 of
# Eq. (8) for the "sequential"/"min" rules and the frozen pre-0.9.9 pipeline
# for the deprecated "legacy" rule.

# Row-wise word proportions as a weighted quanteda dfm and as the plain dense
# matrix the reference implementations consume.
optop_wprop_fixture <- function(fx) {
  wdfm <- quanteda::dfm_weight(quanteda::as.dfm(fx$counts), scheme = "prop")
  W_prop <- sweep(fx$counts, 1, rowSums(fx$counts), "/")
  list(wdfm = wdfm, W_prop = W_prop)
}

# numeric tests run silent: verbose defaults to FALSE and the unconditional
# document-drop alert is muffled where it is not the behavior under test
run_optimal_topic <- function(...) {
  suppressMessages(optimal_topic(..., do_plot = FALSE))
}

expect_matches_reference <- function(got, ref, info = NULL) {
  expect_s3_class(got, "data.table")
  expect_identical(names(got), c("topic", "OpTop", "df", "pval"))
  expect_equal(got$topic, ref$topic, info = info)
  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10, info = info)
  expect_equal(got$df, ref$df, info = info)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10, info = info)
}

test_that("optimal_topic() matches the naive Eq. (8) reference", {
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

test_that("selection = \"legacy\" reproduces the pre-0.9.9 pipeline and warns", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_warning(
    got <- run_optimal_topic(fx$models, wp$wdfm, q = 0.80,
                             selection = "legacy"),
    class = "lifecycle_warning_deprecated"
  )
  ref <- ref_optimal_topic_legacy(fx$models, wp$W_prop, q = 0.80)
  expect_matches_reference(got, ref)
})

test_that("the selection rules pick the expected K", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  ref <- ref_optimal_topic(fx$models, wp$W_prop, q = 0.95)

  # sequential: smallest K the test fails to reject
  seq_pick <- ref$topic[ref$pval > 0.05][1]
  expect_false(is.na(seq_pick))
  expect_message(
    optimal_topic(fx$models, wp$wdfm, do_plot = FALSE, verbose = TRUE),
    sprintf("Optimal model has %d topics", seq_pick)
  )

  # min: global minimum of the standardized statistic
  min_pick <- ref$topic[which.min(ref$OpTop)]
  expect_message(
    optimal_topic(fx$models, wp$wdfm, selection = "min",
                  do_plot = FALSE, verbose = TRUE),
    sprintf("Optimal model has %d topics", min_pick)
  )

  # alpha = 1 makes 'pval > alpha' unsatisfiable: the sequential rule must
  # warn and fall back to the global minimum
  expect_message(
    optimal_topic(fx$models, wp$wdfm, alpha = 1, do_plot = FALSE),
    "falling back to the global minimum"
  )
  expect_message(
    optimal_topic(fx$models, wp$wdfm, alpha = 1, do_plot = FALSE,
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

test_that("optimal_topic() is invariant to document order in the dfm", {
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

test_that("optimal_topic() is silent by default and chatty with verbose = TRUE", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_no_message(optimal_topic(fx$models, wp$wdfm, do_plot = FALSE))
  expect_message(
    optimal_topic(fx$models, wp$wdfm, do_plot = FALSE, verbose = TRUE),
    "Optimal model has"
  )
})

test_that("dropping unmatched documents is signalled even when silent", {
  fx <- optop_test_fixture()

  counts_extra <- rbind(fx$counts, docXX = fx$counts[1, ])
  wdfm_extra <- quanteda::dfm_weight(quanteda::as.dfm(counts_extra),
                                     scheme = "prop")

  expect_message(
    optimal_topic(fx$models, wdfm_extra, do_plot = FALSE),
    "Removed 1 document"
  )
})

test_that("optimal_topic() validates its inputs", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_error(optimal_topic("not a list", wp$wdfm),
               "must be a list")
  expect_error(optimal_topic(fx$models[1L], wp$wdfm),
               "multiple LDA models")
  expect_error(optimal_topic(list(1, 2), wp$wdfm),
               "LDA_VEM")
  expect_error(optimal_topic(fx$models, fx$counts),
               "must be a dfm")
  expect_error(optimal_topic(fx$models, wp$wdfm, q = "a"),
               "q must be a numeric")
  expect_error(optimal_topic(fx$models, wp$wdfm, alpha = "a"),
               "alpha must be a numeric")
  expect_error(optimal_topic(fx$models, wp$wdfm, selection = "typo"),
               "should be one of")
})
