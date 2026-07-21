# The compiled cores receive the dfm as an Armadillo sparse matrix, whose
# constructor requires n_rows * n_cols to fit the index type uword. The
# package compiles with ARMA_64BIT_WORD (src/Makevars), so shapes beyond
# 2^31 - 1 cells must cross the boundary; under 32-bit indices the call
# below dies with "SpMat::init(): requested size is too large". The test is
# cheap by construction: the shape is huge (5e9 cells) but the work is
# O(nnz + J) because only one word column is evaluated.

test_that("the sparse boundary accepts shapes beyond 32-bit indexing", {
  J <- 100000L
  W <- 50000L   # J * W = 5e9 > 2^32 - 1: fails without ARMA_64BIT_WORD
  set.seed(99)
  nnz <- 1000L
  N <- Matrix::sparseMatrix(i = sample.int(J, nnz, replace = TRUE),
                            j = sample.int(W, nnz, replace = TRUE),
                            x = stats::rpois(nnz, 3) + 1,
                            dims = c(J, W))
  N <- methods::as(N, "CsparseMatrix")
  out <- optop_index_word_core(
    matrix(0, 0, 0), N, 0L, 1L,
    as.numeric(Matrix::rowSums(N)), 0.5, 1e-12,
    FALSE, TRUE,
    TRUE, FALSE, FALSE, 1L
  )
  expect_identical(dim(out), c(1L, 6L))
  expect_true(all(is.finite(out)))
})
