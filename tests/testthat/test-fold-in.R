# Direct contracts of the optop_fold_in() methods, engine by engine: each
# returns a document-topic matrix over the new documents with unit row
# sums and the engine's topic count, aligned by document name. The
# holdout layer exercises these indirectly; these tests pin the method
# contracts by name.

test_that("topicmodels fits fold new documents in", {
  hx <- optop_holdout_fixture()
  m <- hx$models[[2]]
  K <- m@k

  th <- optop_fold_in(m, hx$ev)
  expect_true(is.matrix(th))
  expect_identical(dim(th), c(nrow(hx$ev), K))
  expect_identical(rownames(th), rownames(hx$ev))
  expect_equal(unname(rowSums(th)), rep(1, nrow(hx$ev)), tolerance = 1e-8)
  expect_true(all(th >= 0))
})

test_that("seededlda fits fold new documents in", {
  skip_if_not_installed("seededlda")
  hx <- optop_holdout_fixture()
  tr_dfm <- quanteda::as.dfm(as.matrix(hx$tr))
  fit <- seededlda::textmodel_lda(tr_dfm, k = 3, max_iter = 200,
                                  verbose = FALSE)

  th <- optop_fold_in(fit, hx$ev)
  expect_identical(dim(th), c(nrow(hx$ev), 3L))
  expect_identical(rownames(th), rownames(hx$ev))
  expect_equal(unname(rowSums(th)), rep(1, nrow(hx$ev)), tolerance = 1e-8)
})

test_that("WarpLDA wrappers fold new documents in", {
  skip_if_not_installed("text2vec")
  hx <- optop_holdout_fixture()
  lda <- text2vec::LDA$new(n_topics = 3)
  theta_tr <- suppressWarnings(
    lda$fit_transform(hx$tr, n_iter = 100, progressbar = FALSE)
  )
  fit <- optop_warplda(lda, theta_tr)

  th <- suppressWarnings(optop_fold_in(fit, hx$ev, n_iter = 50,
                                       progressbar = FALSE))
  expect_identical(dim(th), c(nrow(hx$ev), 3L))
  expect_identical(rownames(th), rownames(hx$ev))
  expect_equal(unname(rowSums(th)), rep(1, nrow(hx$ev)), tolerance = 1e-8)
})
