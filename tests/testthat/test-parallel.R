# The OpenMP layer of the compiled cores: results must never depend on the
# number of threads (deterministic per-document RNG streams + fixed-order
# reductions), seeds must reproduce, and the compiled bootstrap must target
# the same null distribution as the pure-R rmultinom() oracle.

test_that("the statistic table is bit-identical for any thread count", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  run <- function(threads) {
    suppressMessages(
      optimal_topic(fx$models, wp$wdfm, n_threads = threads,
                    do_plot = FALSE, verbose = FALSE)
    )
  }
  expect_identical(run(4L), run(1L))
})

test_that("bootstrap calibration is bit-identical for any thread count", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  dl <- rowSums(fx$counts)

  run <- function(threads) {
    suppressMessages(
      optimal_topic(fx$models, wp$wdfm, calibrate = "bootstrap",
                    n_boot = 100, doc_lengths = dl, seed = 42,
                    n_threads = threads, do_plot = FALSE, verbose = FALSE)
    )
  }
  expect_identical(run(4L), run(1L))
})

test_that("the compiled bootstrap is seed-reproducible, directly and via set.seed", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  m <- fx$models[[2]]
  env <- optop_envelope(m, wp$wdfm)
  dl <- rowSums(fx$counts)

  # explicit seed
  a <- .optop_boot_null(env$probs, dl, n_boot = 50, seed = 7)
  b <- .optop_boot_null(env$probs, dl, n_boot = 50, seed = 7)
  expect_identical(a, b)

  # NULL seed drawn from the session RNG: set.seed() governs as before
  set.seed(1)
  c1 <- .optop_boot_null(env$probs, dl, n_boot = 50)
  set.seed(1)
  c2 <- .optop_boot_null(env$probs, dl, n_boot = 50)
  expect_identical(c1, c2)

  # different seeds give different draws
  d <- .optop_boot_null(env$probs, dl, n_boot = 50, seed = 8)
  expect_false(identical(a, d))
})

test_that("the compiled bootstrap matches the rmultinom() oracle in distribution", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  m <- fx$models[[2]]
  env <- optop_envelope(m, wp$wdfm)
  dl <- rowSums(fx$counts)

  B <- 5000
  T_cpp <- .optop_boot_null(env$probs, dl, n_boot = B, seed = 99,
                            n_threads = 2L)
  set.seed(2024)
  T_r <- ref_boot_null(env$probs, dl, B)

  # both streams estimate the same null: means within combined MC error,
  # variances within a loose relative band, and the exact Haldane moments
  # bracket both
  se_mean <- sqrt(stats::var(T_cpp) / B + stats::var(T_r) / B)
  expect_lt(abs(mean(T_cpp) - mean(T_r)), 5 * se_mean)
  expect_lt(abs(stats::var(T_cpp) - stats::var(T_r)) / stats::var(T_r), 0.30)

  mm <- .optop_moment_null(env$probs, dl)
  expect_lt(abs(mean(T_cpp) - mm$mu), 5 * stats::sd(T_cpp) / sqrt(B))
  expect_lt(abs(stats::var(T_cpp) - mm$sigma2) / mm$sigma2, 0.15)

  # and the calibrated p-values they imply agree within bootstrap resolution
  T_obs <- env$stat[1L, 2L]
  p_cpp <- (1 + sum(T_cpp >= T_obs)) / (B + 1)
  p_r <- (1 + sum(T_r >= T_obs)) / (B + 1)
  expect_lt(abs(p_cpp - p_r), 5 * sqrt(0.25 / B) + 2 / (B + 1))
})

test_that("n_threads is validated and the OpenMP probe returns a flag", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)

  expect_error(
    optimal_topic(fx$models, wp$wdfm, n_threads = 0,
                  do_plot = FALSE, verbose = FALSE),
    "n_threads"
  )
  expect_error(
    optimal_topic(fx$models, wp$wdfm, n_threads = c(1, 2),
                  do_plot = FALSE, verbose = FALSE),
    "n_threads"
  )

  expect_true(is.logical(optop_openmp_available()))
  expect_length(optop_openmp_available(), 1L)
})

test_that("the index engine and partition kernel are thread-invariant", {
  fx <- optop_test_fixture()

  part4 <- optop_make_partition(fx$models, fx$dtm, c = 5, n_threads = 4L)
  expect_identical(part4$nonrare_offsets, fx$partition$nonrare_offsets)
  expect_identical(part4$nonrare_words, fx$partition$nonrare_words)
  expect_identical(part4$chisq_min_ok, fx$partition$chisq_min_ok)

  tab <- function(threads, lvl) {
    as.data.frame(optop_index_table(fx$models, fx$dtm,
                                    metrics = c("se", "chisq", "deviance"),
                                    macro = TRUE, level = lvl,
                                    partition = fx$partition,
                                    baseline = fx$baseline,
                                    n_threads = threads))
  }
  expect_identical(tab(4L, "document"), tab(1L, "document"))
  expect_identical(tab(4L, "word"), tab(1L, "word"))

  one <- optop_index_deviance(fx$models[[2]], fx$dtm, fx$partition,
                              fx$baseline, macro = TRUE, n_threads = 4L)
  ref <- optop_index_deviance(fx$models[[2]], fx$dtm, fx$partition,
                              fx$baseline, macro = TRUE, n_threads = 1L)
  expect_identical(one, ref)
})
