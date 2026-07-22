# Sharded corpora and lazy models. Every OpTop statistic is a sum of
# independent per-document terms, so a sharded evaluation must reproduce the
# unsharded one: per-document outputs bit-identically, corpus-level sums up
# to floating-point summation order (the grouping of the additions changes
# with the shard boundaries).

corpus_split <- function(m, at) {
  optop_corpus(list(m[1:at, ], m[(at + 1):nrow(m), ]))
}

corpus_index_fun <- function(metric) {
  switch(metric,
         se = optop_index_se,
         chisq = optop_index_chisq,
         deviance = optop_index_deviance)
}

test_that("the constructor validates shards", {
  fx <- optop_test_fixture()
  expect_error(optop_corpus(list()), "at least one shard")
  expect_error(optop_corpus("no/such/file.rds", reader = readRDS),
               "not found")
  expect_error(optop_corpus(character(0)), "at least one shard")
  expect_error(optop_corpus(tempfile()), "reader")
  bad_vocab <- fx$dtm[13:fx$J, rev(seq_len(fx$W))]
  expect_error(optop_corpus(list(fx$dtm[1:12, ], bad_vocab)),
               "same vocabulary")
  expect_error(optop_corpus(list(fx$dtm[1:12, ], fx$dtm[1:12, ])),
               "disjoint")
  corp <- corpus_split(fx$dtm, 12)
  expect_s3_class(corp, "optop_corpus")
  expect_output(print(corp), "2 in-memory shards")
})

test_that("partition, baseline, and indices are shard-invariant", {
  fx <- optop_test_fixture()
  corp <- corpus_split(fx$dtm, 12)

  p1 <- optop_make_partition(fx$models, fx$dtm, c = 5)
  p2 <- optop_make_partition(fx$models, corp, c = 5)
  expect_identical(p1$nonrare_offsets, p2$nonrare_offsets)
  expect_identical(p1$nonrare_words, p2$nonrare_words)
  expect_identical(p1$chisq_min_ok, p2$chisq_min_ok)
  expect_equal(unname(p1$L), unname(p2$L))

  b1 <- optop_make_baseline(fx$dtm)
  b2 <- optop_make_baseline(corp)
  expect_identical(b1$pi_glob, b2$pi_glob)

  # document-level scores are per-document and therefore bit-identical
  for (metric in c("se", "chisq", "deviance")) {
    fun <- corpus_index_fun(metric)
    r1 <- suppressMessages(fun(fx$models[[2]], fx$dtm, p1, b1, macro = TRUE))
    r2 <- suppressMessages(fun(fx$models[[2]], corp, p2, b2, macro = TRUE))
    expect_identical(r1$r2_doc, r2$r2_doc)
    expect_identical(r1$r2, r2$r2)
    # word-level accumulators regroup across shards: summation-order noise
    w1 <- fun(fx$models[[2]], fx$dtm, p1, b1, level = "word")
    w2 <- fun(fx$models[[2]], corp, p2, b2, level = "word")
    expect_equal(w1$r2_word, w2$r2_word, tolerance = 1e-12)
  }

  # the grid table flows through the same engine
  t1 <- suppressMessages(
    optop_index_table(fx$models, fx$dtm, partition = p1, baseline = b1,
                      macro = TRUE)
  )
  t2 <- suppressMessages(
    optop_index_table(fx$models, corp, partition = p2, baseline = b2,
                      macro = TRUE)
  )
  expect_equal(as.data.frame(t1), as.data.frame(t2), tolerance = 1e-12)
})

test_that("lazy model loaders reproduce in-memory grids", {
  fx <- optop_test_fixture()
  loaders <- lapply(fx$models, function(m) {
    force(m)
    function() m
  })
  p1 <- fx$partition
  b1 <- fx$baseline
  p2 <- optop_make_partition(loaders, fx$dtm, c = 5)
  expect_identical(p1$nonrare_offsets, p2$nonrare_offsets)
  expect_identical(p1$nonrare_words, p2$nonrare_words)

  r1 <- suppressMessages(
    optop_index_deviance(fx$models[[2]], fx$dtm, p1, b1, macro = TRUE)
  )
  r2 <- suppressMessages(
    optop_index_deviance(loaders[[2]], fx$dtm, p1, b1, macro = TRUE)
  )
  expect_identical(r1$r2, r2$r2)
  expect_identical(r1$r2_doc, r2$r2_doc)
})

test_that("optop_select is shard-invariant, calibration included", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  loaders <- lapply(fx$models, function(m) {
    force(m)
    function() m
  })
  r_ref <- optop_select(fx$models, wp$wdfm, do_plot = FALSE,
                        verbose = FALSE)
  corp1 <- optop_corpus(wp$wdfm)
  corp2 <- corpus_split(wp$wdfm, 12)

  # the one-shard corpus reproduces the plain path bit for bit
  r1 <- optop_select(fx$models, corp1, do_plot = FALSE, verbose = FALSE)
  expect_identical(as.data.frame(r_ref), as.data.frame(r1))

  r2 <- optop_select(loaders, corp2, do_plot = FALSE, verbose = FALSE)
  expect_equal(r_ref$OpTop, r2$OpTop, tolerance = 1e-12)
  expect_identical(r_ref$df, r2$df)

  dl <- stats::setNames(as.numeric(Matrix::rowSums(fx$dtm)),
                        rownames(fx$counts))
  cb_ref <- optop_select(fx$models, wp$wdfm, calibrate = "bootstrap",
                         n_boot = 40, doc_lengths = dl, seed = 11,
                         do_plot = FALSE, verbose = FALSE)
  cb1 <- optop_select(fx$models, corp1, calibrate = "bootstrap",
                      n_boot = 40, doc_lengths = dl, seed = 11,
                      do_plot = FALSE, verbose = FALSE)
  cb2 <- optop_select(loaders, corp2, calibrate = "bootstrap",
                      n_boot = 40, doc_lengths = dl, seed = 11,
                      do_plot = FALSE, verbose = FALSE)
  expect_identical(cb_ref$pval, cb1$pval)
  expect_equal(cb_ref$pval, cb2$pval, tolerance = 1e-12)

  cm_ref <- optop_select(fx$models, wp$wdfm, calibrate = "moment",
                         doc_lengths = dl, do_plot = FALSE, verbose = FALSE)
  cm2 <- optop_select(loaders, corp2, calibrate = "moment",
                      doc_lengths = dl, do_plot = FALSE, verbose = FALSE)
  expect_equal(cm_ref$OpTop, cm2$OpTop, tolerance = 1e-12)

  # corpus calibration requires named doc_lengths
  expect_error(
    optop_select(fx$models, corp2, calibrate = "moment",
                 doc_lengths = unname(dl), do_plot = FALSE,
                 verbose = FALSE),
    "named doc_lengths"
  )
})

test_that("the bootstrap document streams are keyed by the global index", {
  fx <- optop_test_fixture()
  wp <- optop_wprop_fixture(fx)
  env <- optop_envelope(fx$models[[2]], wp$wdfm, q = 0.95)
  dl <- as.numeric(Matrix::rowSums(fx$dtm))
  J <- fx$J
  t_all <- OpTop:::.optop_boot_null(env$probs, dl, n_boot = 50, seed = 7)
  # every document as its own shard, summed left to right, is the serial
  # reduction of the unsharded core: bit-identical
  per_doc <- lapply(seq_len(J), function(j) {
    OpTop:::.optop_boot_null(env$probs[j], dl[j], n_boot = 50, seed = 7,
                             doc_offset = j - 1)
  })
  expect_identical(Reduce(`+`, per_doc), t_all)
  # arbitrary shard boundaries regroup the additions: ulp-scale noise only
  t_sh <- OpTop:::.optop_boot_null(env$probs[1:12], dl[1:12], n_boot = 50,
                                   seed = 7) +
    OpTop:::.optop_boot_null(env$probs[13:J], dl[13:J], n_boot = 50,
                             seed = 7, doc_offset = 12)
  expect_equal(t_all, t_sh, tolerance = 1e-12)
})

test_that("holdout, gains, and moment tests are shard-invariant", {
  hfx <- optop_holdout_fixture()
  corp <- corpus_split(hfx$ev, 5)
  loaders <- lapply(hfx$models, function(m) {
    force(m)
    function() m
  })

  ho1 <- suppressWarnings(
    optop_index_holdout(hfx$models, hfx$ev, hfx$baseline)
  )
  ho2 <- suppressWarnings(
    optop_index_holdout(loaders, corp, hfx$baseline)
  )
  expect_equal(as.data.frame(ho1$summary), as.data.frame(ho2$summary),
               tolerance = 1e-10)
  for (metric in ho1$metrics) {
    expect_equal(ho1$scores[[metric]], ho2$scores[[metric]],
                 tolerance = 1e-10)
  }

  gt1 <- optop_gain_table(ho1)
  gt2 <- optop_gain_table(ho2)
  expect_equal(as.data.frame(gt1$gains), as.data.frame(gt2$gains),
               tolerance = 1e-10)
  expect_identical(gt1$k_hat, gt2$k_hat)

  mt1 <- suppressWarnings(
    optop_moment_test(hfx$models, hfx$ev, hfx$tr, type = "strata", bins = 4)
  )
  mt2 <- suppressWarnings(
    optop_moment_test(loaders, corp, optop_corpus(list(hfx$tr[1:10, ],
                                                       hfx$tr[11:20, ])),
                      type = "strata", bins = 4)
  )
  expect_equal(as.data.frame(mt1$summary), as.data.frame(mt2$summary),
               tolerance = 1e-10)
})
