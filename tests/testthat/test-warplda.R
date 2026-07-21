# The text2vec WarpLDA adapter: the wrapper carries the user-kept theta and
# the model's topic-word distribution through the standard adapter contract,
# the raw R6 object fails with the documented pointer, and fold-in runs
# through model$transform().

test_that("the WarpLDA wrapper feeds the adapter contract", {
  skip_if_not_installed("text2vec")
  fx <- optop_test_fixture()

  lda <- text2vec::LDA$new(n_topics = 3)
  theta <- suppressWarnings(
    lda$fit_transform(fx$dtm, n_iter = 100, progressbar = FALSE)
  )
  fit <- optop_warplda(lda, theta)

  tp <- optop_as_theta_phi(fit)
  expect_identical(tp$K, 3L)
  expect_identical(tp$docs, rownames(fx$counts))
  expect_identical(tp$terms, colnames(fx$counts))
  expect_equal(unname(rowSums(tp$theta)), rep(1, fx$J), tolerance = 1e-8)
  expect_equal(unname(rowSums(tp$phi)), rep(1, 3), tolerance = 1e-8)

  # the wrapper flows through the partition and the index family
  part <- optop_make_partition(list(fit), fx$dtm, c = 1)
  expect_identical(part$format, 2L)
  expect_length(part$nonrare_offsets, fx$J + 1L)
  res <- optop_index_deviance(fit, fx$dtm, part, fx$baseline, macro = TRUE)
  expect_true(is.finite(res$r2))
  expect_lte(res$r2, 1)
})

test_that("raw WarpLDA objects fail with the wrapper pointer", {
  skip_if_not_installed("text2vec")
  lda <- text2vec::LDA$new(n_topics = 2)
  # the R6 class chain contains "LDA" and "TopicModel": the guard must fire
  # before the topicmodels method can reach posterior()
  expect_error(optop_as_theta_phi(lda), "optop_warplda")
})

test_that("optop_warplda validates its inputs", {
  skip_if_not_installed("text2vec")
  lda <- text2vec::LDA$new(n_topics = 2)
  expect_error(optop_warplda(list(), matrix(1)), "WarpLDA")
  expect_error(optop_warplda(lda, matrix(1, 2, 2)), "row names")
})

test_that("WarpLDA wrappers fold in through the holdout layer", {
  skip_if_not_installed("text2vec")
  hx <- optop_holdout_fixture()

  models <- lapply(2:3, function(k) {
    lda <- text2vec::LDA$new(n_topics = k)
    theta <- suppressWarnings(
      lda$fit_transform(hx$tr, n_iter = 200, progressbar = FALSE)
    )
    optop_warplda(lda, theta)
  })
  ho <- suppressWarnings(
    optop_index_holdout(models, hx$ev, hx$baseline, c = 1,
                        metrics = "deviance",
                        n_iter = 100, progressbar = FALSE)
  )
  expect_identical(dim(ho$scores$deviance), c(10L, 2L))
  expect_true(all(is.finite(ho$summary$macro)))
})

test_that("unbinned aggregation commutes between documents and words", {
  # Lemma 2 of the paper: with the evaluation support equal to the full
  # vocabulary (all-FALSE mask), the total fitted deviance is the same
  # whether summed by document or by word (Poisson-corrected)
  fx <- optop_test_fixture()
  m <- fx$models[[2]]
  unbinned <- list(
    nonrare_offsets = as.numeric(seq(0, fx$J * fx$W, by = fx$W)),
    nonrare_words = rep(seq_len(fx$W) - 1L, fx$J),
    vocab = colnames(fx$counts),
    L = fx$partition$L,
    chisq_min_ok = rep(TRUE, fx$J),
    chisq_min_report = list(n_excluded = 0L, share = 0,
                            excluded_mass = NA_real_),
    c = 1,
    format = 2L
  )
  doc <- optop_index_deviance(m, fx$dtm, unbinned, fx$baseline)
  word <- optop_index_deviance(m, fx$dtm, unbinned, fx$baseline,
                               level = "word")
  expect_equal(sum(doc$d_model), sum(word$d_model), tolerance = 1e-10)
})
