# Alignment tests for the revised methodology (Lewis & Grossetti 2026): the
# harmonized union includes the no-topics baseline, the Pearson min bin obeys
# the grid-wide inclusion rule, aggregation restricts to non-degenerate
# units, and the word-level fitted deviance is the Poisson unit deviance.
# Hand-built nlp_topic_fit objects give exact control over theta and phi
# without fitting.

toy_model <- function(phi_row, docs, vocab) {
  J <- length(docs)
  theta <- matrix(1 / 2, J, 2L, dimnames = list(docs, NULL))
  phi <- rbind(phi_row, phi_row)
  colnames(phi) <- vocab
  structure(list(engine = "toy", model = "lda", method = "toy",
                 model_object = NULL,
                 dtw = theta, tww = phi,
                 doc_ids = docs, vocab = vocab),
            class = c("nlp_topic_fit", "list"))
}

toy_corpus_baseline_rare <- function() {
  # w4 is rare through the baseline only: pi_glob(w4) = 1/60 < tau = 0.05,
  # while every fitted probability stays >= 0.2 for both models
  vocab <- c("w1", "w2", "w3", "w4")
  docs <- c("d1", "d2", "d3")
  counts <- rbind(c(7L, 6L, 6L, 1L),
                  c(7L, 7L, 6L, 0L),
                  c(7L, 7L, 6L, 0L))
  dimnames(counts) <- list(docs, vocab)
  dtm <- methods::as(counts, "CsparseMatrix")
  models <- list(toy_model(c(0.25, 0.25, 0.25, 0.25), docs, vocab),
                 toy_model(c(0.30, 0.30, 0.20, 0.20), docs, vocab))
  list(dtm = dtm, counts = counts, models = models,
       docs = docs, vocab = vocab)
}

test_that("the no-topics baseline enters the harmonized union", {
  toy <- toy_corpus_baseline_rare()
  part <- optop_make_partition(toy$models, toy$dtm, c = 1)

  # every fitted probability clears tau = 0.05, so without the baseline term
  # nothing would be rare; pi_glob(w4) = 1/60 makes w4 rare in every document
  expect_true(all(part$rare_mask[, "w4"]))
  expect_false(any(part$rare_mask[, c("w1", "w2", "w3")]))
})

test_that("the Pearson min bin is excluded when the baseline mass is small", {
  toy <- toy_corpus_baseline_rare()
  part <- optop_make_partition(toy$models, toy$dtm, c = 1)

  # B_min = L_j * pi_glob(w4) = 20/60 < 1 while min_K E_min = 20 * 0.2 >= 1:
  # the rule excludes the min bin for every document
  expect_false(any(part$chisq_min_ok))
  expect_identical(part$chisq_min_report$n_excluded, 3L)
  expect_equal(part$chisq_min_report$share, 1)
  expect_equal(part$chisq_min_report$excluded_mass, mean(c(1 / 20, 0, 0)))

  baseline <- optop_make_baseline(toy$dtm)
  m <- toy$models[[1]]
  expect_message(
    got <- optop_index_chisq(m, toy$dtm, part, baseline, macro = TRUE),
    "Pearson min bin excluded"
  )

  # naive value with the min bin dropped from BOTH sides, per document
  tp <- optop_as_theta_phi(m)
  pi_row <- baseline$pi_glob[toy$vocab]
  keep <- c("w1", "w2", "w3")
  for (j in 1:3) {
    L_j <- sum(toy$counts[j, ])
    E_j <- L_j * drop(tp$theta[j, ] %*% tp$phi)[keep]
    B_j <- L_j * pi_row[keep]
    N_j <- toy$counts[j, keep]
    dK <- sum((N_j - E_j)^2 / E_j)
    dnull <- sum((N_j - B_j)^2 / B_j)
    expect_equal(got$d_model[j], unname(dK), tolerance = 1e-12)
    expect_equal(got$d_null[j], unname(dnull), tolerance = 1e-12)
  }

  # deviance and squared error keep the min bin: their per-document nulls
  # must include the w4 term, hence exceed the keep-only computation for the
  # document that observed w4
  dv <- optop_index_deviance(m, toy$dtm, part, baseline)
  expect_gt(dv$d_null[1],
            2 * sum(toy$counts[1, keep] *
                      log(toy$counts[1, keep] /
                            (sum(toy$counts[1, ]) * pi_row[keep]))))
})

test_that("the inclusion rule keeps a min bin with enough grid-wide mass", {
  # two model-rare words with fitted mass 0.008 each (tau = 0.01) and
  # baseline mass 0.02 each: E_min = 100 * 0.016 and B_min = 100 * 0.04 both
  # clear c = 1, so the min bin stays in
  vocab <- c("w1", "w2", "w3", "w4", "w5")
  docs <- c("d1", "d2", "d3")
  counts <- rbind(c(32L, 31L, 33L, 2L, 2L),
                  c(31L, 32L, 33L, 2L, 2L),
                  c(31L, 31L, 34L, 2L, 2L))
  dimnames(counts) <- list(docs, vocab)
  dtm <- methods::as(counts, "CsparseMatrix")
  phi_row <- c(0.40, 0.30, 0.284, 0.008, 0.008)
  models <- list(toy_model(phi_row, docs, vocab))

  part <- optop_make_partition(models, dtm, c = 1)
  expect_true(all(part$rare_mask[, c("w4", "w5")]))
  expect_false(any(part$rare_mask[, c("w1", "w2", "w3")]))
  expect_true(all(part$chisq_min_ok))
  expect_identical(part$chisq_min_report$n_excluded, 0L)

  # no exclusion message when every min bin passes the rule; the three
  # near-identical documents all sit close to the pooled baseline, so the
  # 0.14.1 null-discrepancy floor is disabled to isolate the min-bin path
  baseline <- optop_make_baseline(dtm)
  expect_no_message(
    optop_index_chisq(models[[1]], dtm, part, baseline, min_null = 0)
  )
})

test_that("degenerate documents are excluded from the Micro index", {
  # doc d3 reproduces the supplied baseline exactly, so its null discrepancy
  # is zero for every family while its fitted discrepancy is not
  vocab <- c("w1", "w2", "w3", "w4")
  docs <- c("d1", "d2", "d3")
  counts <- rbind(c(8L, 5L, 4L, 3L),
                  c(2L, 7L, 6L, 5L),
                  c(2L, 3L, 4L, 1L))
  dimnames(counts) <- list(docs, vocab)
  dtm <- methods::as(counts, "CsparseMatrix")
  models <- list(toy_model(c(0.25, 0.25, 0.25, 0.25), docs, vocab),
                 toy_model(c(0.30, 0.30, 0.20, 0.20), docs, vocab))
  part <- optop_make_partition(models, dtm, c = 1)

  # user-supplied baseline equal to d3's empirical distribution
  baseline <- list(pi_glob = stats::setNames(c(0.2, 0.3, 0.4, 0.1), vocab))

  for (fun in list(optop_index_se, optop_index_chisq, optop_index_deviance)) {
    got <- fun(models[[1]], dtm, part, baseline, macro = TRUE)
    expect_identical(got$d_null[3], 0)
    expect_gt(got$d_model[3], 0)
    expect_true(is.na(got$r2_doc[3]))
    # Micro over J+ only: the degenerate document's fitted discrepancy must
    # not leak into the numerator
    expect_equal(got$r2,
                 1 - sum(got$d_model[1:2]) / sum(got$d_null[1:2]),
                 tolerance = 1e-12)
    expect_false(isTRUE(all.equal(
      got$r2, 1 - sum(got$d_model) / sum(got$d_null)
    )))
    expect_equal(got$r2_macro, mean(got$r2_doc[1:2]), tolerance = 1e-12)
  }
})

test_that("the word-level fitted deviance is the Poisson unit deviance", {
  fx <- optop_test_fixture()
  m <- fx$models[[2]]
  got <- optop_index_deviance(m, fx$dtm, fx$partition, fx$baseline,
                              level = "word")

  # unit deviances are nonnegative, so the index never exceeds one
  expect_true(all(got$d_model >= -1e-9))
  expect_true(all(got$r2_word <= 1 + 1e-12, na.rm = TRUE))

  # the correction changes the fitted side whenever counts and expectations
  # disagree in total, which the uncorrected formula ignored
  tp <- ref_theta_phi(m)
  E <- ref_expected(tp$theta, tp$phi, fx$partition$L)
  N <- as.matrix(fx$dtm)
  uncorrected <- sapply(seq_len(fx$W), function(w) {
    idx <- N[, w] > 0
    if (!any(idx)) return(0)
    2 * sum(N[idx, w] * (log(N[idx, w]) - log(pmax(E[idx, w], 1e-12))))
  })
  correction <- 2 * (colSums(N) - colSums(E))
  expect_equal(unname(got$d_model), unname(uncorrected - correction),
               tolerance = 1e-10)
  expect_gt(max(abs(correction)), 0)
})

test_that("partitions from OpTop < 0.13.0 fall back with a warning", {
  fx <- optop_test_fixture()
  m <- fx$models[[1]]
  old_part <- fx$partition[c("rare_mask", "L")]

  expect_warning(
    old <- optop_index_chisq(m, fx$dtm, old_part, fx$baseline),
    "chisq_min_ok"
  )

  # the fallback keeps every min bin, i.e. the all-TRUE behavior
  forced <- fx$partition
  forced$chisq_min_ok <- rep(TRUE, fx$J)
  forced$chisq_min_report <- list(n_excluded = 0L, share = 0,
                                  excluded_mass = NA_real_)
  new <- optop_index_chisq(m, fx$dtm, forced, fx$baseline)
  expect_equal(old$r2, new$r2, tolerance = 1e-15)
  expect_equal(old$d_model, new$d_model, tolerance = 1e-15)

  # deviance and SE need no flag and stay silent
  expect_no_warning(optop_index_deviance(m, fx$dtm, old_part, fx$baseline))
  expect_no_warning(optop_index_se(m, fx$dtm, old_part, fx$baseline))
})

test_that("the non-paper features are deprecated", {
  fx <- optop_test_fixture()
  m <- fx$models[[1]]

  lifecycle::expect_deprecated(
    optop_index_se(m, fx$dtm, fx$partition, fx$baseline,
                   macro = TRUE, ztest = TRUE),
    "ztest"
  )
  lifecycle::expect_deprecated(
    optop_index_se(m, fx$dtm, fx$partition, fx$baseline, reopt = "se"),
    "reopt"
  )
  lifecycle::expect_deprecated(
    optop_index_se(m, fx$dtm, fx$partition, fx$baseline,
                   add_baseline_topic = TRUE),
    "add_baseline_topic"
  )
  lifecycle::expect_deprecated(
    optop_index_chisq(m, fx$dtm, fx$partition, fx$baseline,
                      macro = TRUE, ztest = TRUE),
    "ztest"
  )
  lifecycle::expect_deprecated(
    optop_index_deviance(m, fx$dtm, fx$partition, fx$baseline,
                         macro = TRUE, ztest = TRUE),
    "ztest"
  )
  lifecycle::expect_deprecated(
    optop_index_table(fx$models[1:2], fx$dtm, metrics = "se",
                      partition = fx$partition, baseline = fx$baseline,
                      ztest = TRUE),
    "ztest"
  )
})

test_that("default c is 1 in partition and table", {
  expect_identical(formals(optop_make_partition)$c, 1)
  expect_identical(formals(optop_index_table)$c, 1)
  expect_identical(formals(optop_index_se)$add_baseline_topic, FALSE)
  expect_identical(formals(optop_index_chisq)$add_baseline_topic, FALSE)
  expect_identical(formals(optop_index_table)$add_baseline_topic, FALSE)
})
