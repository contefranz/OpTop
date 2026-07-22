# The self-contained LDA corpus simulator: exact document lengths, seed
# reproducibility, dimnames, the Dirichlet path, the input gates, and a
# large-sample agreement with the mixture probabilities it draws from.

sim_tp <- function(J = 12L, k = 3L, W = 40L, seed = 5) {
  set.seed(seed)
  rdirich <- function(n, m) {
    g <- matrix(stats::rgamma(n * m, shape = 1), n, m)
    g / rowSums(g)
  }
  theta <- rdirich(J, k)
  rownames(theta) <- sprintf("doc%02d", seq_len(J))
  phi <- rdirich(k, W)
  colnames(phi) <- sprintf("w%02d", seq_len(W))
  list(theta = theta, phi = phi, J = J, k = k, W = W)
}

test_that("sim_dfm draws a counts dfm with exact document lengths", {
  tp <- sim_tp()
  L <- sample(50:150, tp$J, replace = TRUE)
  out <- sim_dfm(tp$theta, tp$phi, doc_length = L, seed = 1)

  expect_s4_class(out, "dfm")
  expect_identical(dim(out), c(tp$J, tp$W))
  expect_identical(quanteda::docnames(out), rownames(tp$theta))
  expect_identical(quanteda::featnames(out), colnames(tp$phi))
  expect_equal(unname(Matrix::rowSums(out)), as.numeric(L))
  expect_true(all(as.matrix(out) >= 0))
  expect_true(all(as.matrix(out) == round(as.matrix(out))))
})

test_that("a scalar doc_length is recycled and default names are generated", {
  tp <- sim_tp()
  theta <- unname(tp$theta)
  phi <- unname(tp$phi)
  out <- sim_dfm(theta, phi, doc_length = 80, seed = 2)

  expect_identical(dim(out), c(tp$J, tp$W))
  expect_equal(unname(Matrix::rowSums(out)), rep(80, tp$J))
  expect_identical(quanteda::docnames(out), sprintf("text%d", seq_len(tp$J)))
  expect_identical(quanteda::featnames(out), sprintf("feat%d", seq_len(tp$W)))
})

test_that("sim_dfm is seed-reproducible", {
  tp <- sim_tp()
  a <- sim_dfm(tp$theta, tp$phi, doc_length = 100, seed = 7)
  b <- sim_dfm(tp$theta, tp$phi, doc_length = 100, seed = 7)
  d <- sim_dfm(tp$theta, tp$phi, doc_length = 100, seed = 8)
  expect_identical(as.matrix(a), as.matrix(b))
  expect_false(identical(as.matrix(a), as.matrix(d)))

  # NULL seed draws from the session RNG, governed by set.seed() as usual
  set.seed(3)
  c1 <- sim_dfm(tp$theta, tp$phi, doc_length = 100)
  set.seed(3)
  c2 <- sim_dfm(tp$theta, tp$phi, doc_length = 100)
  expect_identical(as.matrix(c1), as.matrix(c2))
})

test_that("the alpha path draws fresh Dirichlet topic weights", {
  tp <- sim_tp()
  out <- sim_dfm(tp$theta, tp$phi, doc_length = 120, alpha = 0.5, seed = 4)
  expect_identical(dim(out), c(tp$J, tp$W))
  expect_equal(unname(Matrix::rowSums(out)), rep(120, tp$J))
  expect_identical(quanteda::docnames(out), rownames(tp$theta))

  # a vector alpha of length k is accepted; other lengths are not
  out_vec <- sim_dfm(tp$theta, tp$phi, doc_length = 120,
                     alpha = rep(1, tp$k), seed = 4)
  expect_identical(dim(out_vec), c(tp$J, tp$W))
  expect_error(
    sim_dfm(tp$theta, tp$phi, doc_length = 120, alpha = c(1, 2)),
    "alpha"
  )

  # with alpha the supplied DTW rows need not be distributions
  theta_raw <- matrix(5, tp$J, tp$k)
  out_raw <- sim_dfm(theta_raw, tp$phi, doc_length = 50, alpha = 1, seed = 9)
  expect_equal(unname(Matrix::rowSums(out_raw)), rep(50, tp$J))
})

test_that("large samples agree with the mixture probabilities", {
  tp <- sim_tp(J = 6L)
  L <- 50000
  out <- sim_dfm(tp$theta, tp$phi, doc_length = L, seed = 10)
  P <- tp$theta %*% tp$phi
  prop <- as.matrix(out) / L
  # multinomial noise at L = 5e4 keeps every cell within a loose band
  expect_lt(max(abs(prop - P)), 6 * sqrt(max(P) / L))
})

test_that("sim_dfm validates its inputs", {
  tp <- sim_tp()
  expect_error(sim_dfm("x", tp$phi, 10), "matrix or a data.frame")
  expect_error(sim_dfm(tp$theta, "x", 10), "matrix or a data.frame")
  expect_error(sim_dfm(tp$theta, tp$phi, -5), "positive")
  expect_error(sim_dfm(tp$theta, tp$phi, c(10, 20)),
               "one entry per document")
  expect_error(sim_dfm(tp$theta, tp$phi[-1, ], 10), "one row per topic")
  expect_error(sim_dfm(tp$theta * 2, tp$phi, 10), "summing to 1")
  expect_error(sim_dfm(tp$theta, tp$phi * 2, 10), "summing to 1")

  # data.frame inputs are coerced
  out <- sim_dfm(as.data.frame(tp$theta), as.data.frame(tp$phi),
                 doc_length = 30, seed = 11)
  expect_identical(dim(out), c(tp$J, tp$W))
})
