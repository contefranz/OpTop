# Characterization tests for optimal_topic(): the C++ core must agree with the
# naive per-document reference implementation (helper-reference.R) that mirrors
# the current semantics verbatim, including the round(cum * 1e4) truncation and
# the lower-tail chi-square p-value documented in AUDIT.md.

# Row-wise word proportions as a weighted quanteda dfm and as the plain dense
# matrix the reference implementation consumes.
optop_wprop_fixture <- function(fx) {
  wdfm <- quanteda::dfm_weight(quanteda::as.dfm(fx$counts), scheme = "prop")
  W_prop <- sweep(fx$counts, 1, rowSums(fx$counts), "/")
  list(wdfm = wdfm, W_prop = W_prop)
}

run_optimal_topic <- function(...) {
  res <- NULL
  capture.output(res <- optimal_topic(..., do_plot = FALSE))
  res
}

test_that("optimal_topic() matches the naive reference", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  for (q in c(0.80, 0.95)) {
    got <- run_optimal_topic(fx$models, wp$wdfm, q = q)
    ref <- ref_optimal_topic(fx$models, wp$W_prop, q = q)

    info <- sprintf("q=%.2f", q)
    expect_s3_class(got, "data.table")
    expect_identical(names(got), c("topic", "OpTop", "pval"))
    expect_equal(got$topic, ref$topic, info = info)
    expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10, info = info)
    expect_equal(got$pval, ref$pval, tolerance = 1e-10, info = info)
  }
})

test_that("documents missing from the models are dropped, order preserved", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  counts_extra <- rbind(fx$counts, docXX = fx$counts[1, ])
  wdfm_extra <- quanteda::dfm_weight(quanteda::as.dfm(counts_extra),
                                     scheme = "prop")

  got <- run_optimal_topic(fx$models, wdfm_extra)
  ref <- ref_optimal_topic(fx$models, wp$W_prop)

  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10)
})

test_that("optimal_topic() is invariant to document order in the dfm", {
  skip("document-order alignment fix pending (see AUDIT.md)")

  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  set.seed(42)
  perm <- sample(nrow(fx$counts))
  wdfm_perm <- quanteda::dfm_weight(quanteda::as.dfm(fx$counts[perm, ]),
                                    scheme = "prop")

  got <- run_optimal_topic(fx$models, wdfm_perm)
  ref <- ref_optimal_topic(fx$models, wp$W_prop)

  expect_equal(got$OpTop, ref$OpTop, tolerance = 1e-10)
  expect_equal(got$pval, ref$pval, tolerance = 1e-10)
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
})
