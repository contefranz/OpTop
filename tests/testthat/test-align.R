# The exported alignment helper: the remediation path the index functions
# point at when the dtm and the models disagree on the vocabulary.

test_that("optop_align_dtm_to_models reorders and subsets to the model vocabulary", {
  fx <- optop_test_fixture()
  tp <- ref_theta_phi(fx$models[[2]])
  m <- optop_model(tp$theta, tp$phi)

  # shuffle the columns and add a feature no model knows
  set.seed(1)
  perm <- sample.int(fx$W)
  shuffled <- cbind(fx$counts[, perm], zzz_extra = 1L)
  aligned <- optop_align_dtm_to_models(shuffled, list(m))

  expect_identical(colnames(aligned), m$terms)
  expect_equal(unname(as.matrix(aligned)), unname(fx$counts))

  # the aligned matrix flows through partition, baseline, and the index
  dtm <- methods::as(Matrix::Matrix(as.matrix(aligned), sparse = TRUE),
                     "CsparseMatrix")
  part <- optop_make_partition(list(m), dtm, c = 5)
  base <- optop_make_baseline(dtm)
  res <- suppressMessages(optop_index_deviance(m, dtm, part, base))
  expect_true(is.finite(res$r2))

  # while the shuffled matrix is refused with the pointer here
  shuf_sparse <- methods::as(Matrix::Matrix(shuffled, sparse = TRUE),
                             "CsparseMatrix")
  expect_error(
    optop_index_deviance(m, shuf_sparse, part, base),
    "optop_align_dtm_to_models"
  )
})

test_that("the intersection follows the first model's term order", {
  fx <- optop_test_fixture()
  tp <- ref_theta_phi(fx$models[[2]])

  m_full <- optop_model(tp$theta, tp$phi)
  # a second model missing the first five terms
  keep <- seq(6, fx$W)
  phi_sub <- tp$phi[, keep, drop = FALSE]
  phi_sub <- phi_sub / rowSums(phi_sub)
  m_sub <- optop_model(tp$theta, phi_sub)

  aligned <- optop_align_dtm_to_models(fx$counts, list(m_full, m_sub))
  expect_identical(colnames(aligned), m_full$terms[keep])

  # dfm input returns a dgCMatrix on the common vocabulary
  dfm_in <- quanteda::as.dfm(fx$counts)
  aligned_dfm <- optop_align_dtm_to_models(dfm_in, list(m_full, m_sub))
  expect_s4_class(aligned_dfm, "dgCMatrix")
  expect_identical(colnames(aligned_dfm), m_full$terms[keep])

  # no overlap at all is an error
  phi_alien <- tp$phi
  colnames(phi_alien) <- paste0("alien_", seq_len(fx$W))
  m_alien <- optop_model(tp$theta, phi_alien)
  expect_error(optop_align_dtm_to_models(fx$counts, list(m_full, m_alien)),
               "Empty vocabulary intersection")

  # a plain matrix missing common terms is an error
  expect_error(
    optop_align_dtm_to_models(fx$counts[, 1:10], list(m_full)),
    "missing terms"
  )
})
