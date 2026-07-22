# The public model constructor: validated construction, identity through
# the adapter generic, full flow through partition and indices, and the
# documented no-engine limitation (no fold-in, so no held-out path).

test_that("optop_model builds the adapter contract from bare matrices", {
  fx <- optop_test_fixture()
  tp_ref <- ref_theta_phi(fx$models[[2]])

  m <- optop_model(tp_ref$theta, tp_ref$phi)
  expect_s3_class(m, "optop_theta_phi")
  expect_identical(m$K, ncol(tp_ref$theta))
  expect_identical(m$docs, rownames(tp_ref$theta))
  expect_identical(m$terms, colnames(tp_ref$phi))

  # the adapter generic returns the contract unchanged
  tp <- optop_as_theta_phi(m)
  expect_identical(tp$theta, tp_ref$theta)
  expect_identical(tp$phi, tp_ref$phi)
})

test_that("optop_model flows through the partition and the index family", {
  fx <- optop_test_fixture()
  tps <- lapply(fx$models, ref_theta_phi)
  models <- lapply(tps, function(tp) optop_model(tp$theta, tp$phi))

  part <- optop_make_partition(models, fx$dtm, c = 5)
  expect_identical(part$nonrare_offsets, fx$partition$nonrare_offsets)
  expect_identical(part$nonrare_words, fx$partition$nonrare_words)

  res_model <- suppressMessages(
    optop_index_deviance(models[[2]], fx$dtm, fx$partition, fx$baseline,
                         macro = TRUE)
  )
  res_fit <- suppressMessages(
    optop_index_deviance(fx$models[[2]], fx$dtm, fx$partition, fx$baseline,
                         macro = TRUE)
  )
  expect_identical(res_model$r2, res_fit$r2)
  expect_identical(res_model$r2_doc, res_fit$r2_doc)
})

test_that("optop_model validates the contract on construction", {
  fx <- optop_test_fixture()
  tp <- ref_theta_phi(fx$models[[2]])

  # topic dimensions must agree
  expect_error(optop_model(tp$theta, tp$phi[-1, , drop = FALSE]),
               "number of topics")
  # document identifiers are required, never positional
  theta_anon <- tp$theta
  rownames(theta_anon) <- NULL
  expect_error(optop_model(theta_anon, tp$phi), "document identifier")
  # terms are required
  phi_anon <- tp$phi
  colnames(phi_anon) <- NULL
  expect_error(optop_model(tp$theta, phi_anon), "one term per column")
  # probabilities must be nonnegative and rows must sum to 1
  theta_neg <- tp$theta
  theta_neg[1, 1] <- -theta_neg[1, 1]
  expect_error(optop_model(theta_neg, tp$phi), "negative")
  expect_error(optop_model(2 * tp$theta, tp$phi), "sum to 1")
})

test_that("optop_model prints its header and cannot fold in", {
  fx <- optop_test_fixture()
  tp <- ref_theta_phi(fx$models[[2]])
  m <- optop_model(tp$theta, tp$phi)

  expect_snapshot(print(m))

  # no engine behind the object: the held-out fold-in must refuse it with
  # the pointed message, not the generic default
  expect_error(optop_fold_in(m, fx$dtm), "bare optop_model")
})
