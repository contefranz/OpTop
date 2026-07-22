# The calibration layer: the C++ envelope export against a naive reference,
# the bootstrap null as the oracle for the closed-form moments, and the
# headline property — data simulated from a fitted model gets extreme
# asymptotic chi-square p-values but interior calibrated ones.

# counts drawn from the fitted model itself: the conditional null is true by
# construction, with the fixture's own document names and vocabulary
optop_self_generated <- function(model, counts, seed) {
  I <- model@gamma %*% exp(model@beta)
  set.seed(seed)
  sim <- t(vapply(seq_len(nrow(counts)), function(j) {
    as.integer(stats::rmultinom(1, size = sum(counts[j, ]), prob = I[j, ]))
  }, integer(ncol(counts))))
  dimnames(sim) <- dimnames(counts)
  quanteda::dfm_weight(quanteda::as.dfm(sim), scheme = "prop")
}

test_that("the exported envelope satisfies its invariants and matches the reference", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  for (m in fx$models) {
    env <- optop_envelope(m, wp$wdfm)

    # bin probabilities are a distribution per document
    expect_true(all(abs(vapply(env$probs, sum, numeric(1)) - 1) < 1e-8))

    # sum_j (B_j - 1) is the df the statistic core reports
    expect_equal(sum(lengths(env$probs) - 1L), env$stat[1, 3])

    # the naive R envelope agrees bin by bin
    ref <- ref_envelope(m, q = 0.95)
    expect_equal(env$probs, ref, tolerance = 1e-10,
                 ignore_attr = TRUE)
  }
})

test_that("self-generated data: asymptotic p-values saturate, calibrated ones are interior", {
  fx <- optop_test_fixture()
  gen <- fx$models[[2]]  # k = 3 generates the data
  wdfm_null <- optop_self_generated(gen, fx$counts, seed = 7)
  dl <- rowSums(fx$counts)

  got <- suppressMessages(
    optop_select(fx$models, wdfm_null, calibrate = "bootstrap",
                 n_boot = 500, doc_lengths = dl, seed = 11,
                 do_plot = FALSE, verbose = FALSE)
  )

  expect_identical(names(got),
                   c("topic", "OpTop", "df", "pval", "pval_chisq"))

  row_gen <- got[got$topic == gen@k]
  # the chi-square yardstick is far off scale even under the true model
  expect_true(row_gen$pval_chisq > 0.999 || row_gen$pval_chisq < 0.001)
  # the calibrated p-value lands in the interior; the bounds are deliberately
  # loose because the fixture fits (and hence the null draw) vary slightly
  # across platforms and BLAS builds
  expect_gt(row_gen$pval, 0.001)
  expect_lt(row_gen$pval, 0.999)
})

test_that("bootstrap calibration is reproducible under a seed", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  dl <- rowSums(fx$counts)

  run <- function() {
    suppressMessages(
      optop_select(fx$models, wp$wdfm, calibrate = "bootstrap",
                   n_boot = 100, doc_lengths = dl, seed = 42,
                   do_plot = FALSE, verbose = FALSE)
    )
  }
  expect_identical(run(), run())
})

test_that("the closed-form moments match the bootstrap oracle", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  m <- fx$models[[2]]
  env <- optop_envelope(m, wp$wdfm)
  dl <- rowSums(fx$counts)

  mm <- .optop_moment_null(env$probs, dl)

  set.seed(123)
  B <- 5000
  T_null <- .optop_boot_null(env$probs, dl, B)

  # mean within 5 standard errors, variance within 15% relative error
  expect_lt(abs(mean(T_null) - mm$mu), 5 * stats::sd(T_null) / sqrt(B))
  expect_lt(abs(stats::var(T_null) - mm$sigma2) / mm$sigma2, 0.15)
})

test_that("moment calibration tracks the bootstrap on self-generated data", {
  fx <- optop_test_fixture()
  gen <- fx$models[[2]]
  wdfm_null <- optop_self_generated(gen, fx$counts, seed = 7)
  dl <- rowSums(fx$counts)

  got_mm <- suppressMessages(
    optop_select(fx$models, wdfm_null, calibrate = "moment",
                 doc_lengths = dl, do_plot = FALSE, verbose = FALSE)
  )
  row_gen <- got_mm[got_mm$topic == gen@k]
  expect_gt(row_gen$pval, 0.001)
  expect_lt(row_gen$pval, 0.999)
})

test_that("calibration validates its inputs", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  dl <- rowSums(fx$counts)

  expect_error(
    optop_select(fx$models, wp$wdfm, calibrate = "bootstrap",
                 do_plot = FALSE, verbose = FALSE),
    "needs doc_lengths"
  )
  expect_error(
    optop_select(fx$models, wp$wdfm, calibrate = "bootstrap",
                 doc_lengths = unname(dl)[-1], do_plot = FALSE,
                 verbose = FALSE),
    "one entry per document"
  )
  expect_error(
    optop_select(fx$models, wp$wdfm, calibrate = "bootstrap",
                 doc_lengths = -dl, do_plot = FALSE, verbose = FALSE),
    "positive"
  )
  dl_bad <- stats::setNames(dl, paste0("x", seq_along(dl)))
  expect_error(
    optop_select(fx$models, wp$wdfm, calibrate = "bootstrap",
                 doc_lengths = dl_bad, do_plot = FALSE, verbose = FALSE),
    "missing entries"
  )
  expect_error(
    optop_select(fx$models, wp$wdfm, calibrate = "bootstrap",
                 doc_lengths = dl, n_boot = 5, do_plot = FALSE,
                 verbose = FALSE),
    "n_boot"
  )
})
