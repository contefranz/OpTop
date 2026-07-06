# Branch coverage of the compiled bootstrap sampler (optop_boot_null_core):
# both sampling regimes and their guards, validated against the exact Haldane
# moments of the conditional multinomial null. The draws are deterministic
# given the seed, and the moment tolerances are wide enough to hold for any
# standard-library RNG stream, so none of these tests can flake.

# a normalized probability vector mixing masses above and below 1/k, so the
# alias-table construction exercises both the small and the large stacks;
# bounded away from zero to keep the null moments well conditioned
mixed_probs <- function(k, seed) {
  set.seed(seed)
  p <- c(rgamma(floor(k / 2), shape = 5), rgamma(ceiling(k / 2), shape = 0.5))
  p <- pmax(p, 0.2 * mean(p))
  p / sum(p)
}

test_that("the wide-envelope alias sampler targets the exact null moments", {
  k <- 40L
  p <- mixed_probs(k, seed = 1)
  N <- 12
  B <- 4000L

  # k > N selects the alias-table path
  T_null <- optop_boot_null_core(p, k, N, B, seed = 11, n_threads = 1L)
  expect_length(T_null, B)
  expect_true(all(is.finite(T_null) & T_null >= 0))

  mm <- .optop_moment_null(list(p), N)
  expect_lt(abs(mean(T_null) - mm$mu), 5 * stats::sd(T_null) / sqrt(B))
  expect_lt(abs(stats::var(T_null) - mm$sigma2) / mm$sigma2, 0.30)

  # the thread contract holds on this path too
  expect_identical(
    optop_boot_null_core(p, k, N, B, seed = 11, n_threads = 3L),
    T_null
  )
})

test_that("wide and narrow regimes combine in one statistic", {
  p_wide <- mixed_probs(50L, seed = 2)
  p_narrow <- mixed_probs(5L, seed = 3)
  probs <- c(p_wide, p_narrow)
  counts <- c(50L, 5L)
  lengths <- c(10, 200)
  B <- 4000L

  T_null <- optop_boot_null_core(probs, counts, lengths, B,
                                 seed = 21, n_threads = 2L)
  expect_length(T_null, B)

  mm <- .optop_moment_null(list(p_wide, p_narrow), lengths)
  expect_lt(abs(mean(T_null) - mm$mu), 5 * stats::sd(T_null) / sqrt(B))
  expect_lt(abs(stats::var(T_null) - mm$sigma2) / mm$sigma2, 0.30)
})

test_that("the narrow-path early exit reproduces and keeps the null mean", {
  # a dominant first bin exhausts the count on most replicates, so the
  # remaining bins take the closed-form suffix branch
  k <- 31L
  p <- c(0.99, rep(0.01 / 30, 30))
  N <- 5
  B <- 4000L

  T_null <- optop_boot_null_core(p, k, N, B, seed = 31, n_threads = 1L)
  expect_true(all(is.finite(T_null) & T_null >= 0))

  # the mean check is robust to the heavy tail: the standard error is
  # estimated from the same sample, so a wide tail widens the band
  mm <- .optop_moment_null(list(p), N)
  expect_lt(abs(mean(T_null) - mm$mu), 5 * stats::sd(T_null) / sqrt(B))

  expect_identical(
    optop_boot_null_core(p, k, N, B, seed = 31, n_threads = 1L),
    T_null
  )
})

test_that("the core validates its inputs", {
  p <- mixed_probs(40L, seed = 4)

  expect_error(
    optop_boot_null_core(p, 40L, c(100, 50), 10L, 1, 1L),
    "doc_lengths"
  )
  expect_error(
    optop_boot_null_core(p, 40L, 100, 0L, 1, 1L),
    "n_boot"
  )
  expect_error(
    optop_boot_null_core(p, c(40L, 0L), c(100, 50), 10L, 1, 1L),
    "envelope bin"
  )
  expect_error(
    optop_boot_null_core(p, c(40L, 2L), c(100, 50), 10L, 1, 1L),
    "disagree"
  )
})

test_that("a non-positive thread count is clamped to one thread", {
  p <- mixed_probs(40L, seed = 5)
  expect_identical(
    optop_boot_null_core(p, 40L, 12, 50L, 3, 0L),
    optop_boot_null_core(p, 40L, 12, 50L, 3, 1L)
  )
})

test_that("corpora larger than one document block reduce in fixed order", {
  n_docs <- 1100L  # two blocks of 1024
  probs <- rep(c(0.7, 0.3), n_docs)
  counts <- rep(2L, n_docs)
  lengths <- rep(4, n_docs)

  T_1 <- optop_boot_null_core(probs, counts, lengths, 5L,
                              seed = 41, n_threads = 1L)
  expect_length(T_1, 5L)
  expect_true(all(is.finite(T_1) & T_1 >= 0))
  expect_identical(
    optop_boot_null_core(probs, counts, lengths, 5L,
                         seed = 41, n_threads = 4L),
    T_1
  )
})
