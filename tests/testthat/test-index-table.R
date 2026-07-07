test_that("optop_index_table reproduces the single-model document indices", {
  withr::local_options(lifecycle_verbosity = "quiet")
  fx <- optop_test_fixture()
  tab <- optop_index_table(fx$models, fx$dtm,
                           metrics = c("se", "chisq", "deviance"),
                           macro = TRUE, ztest = TRUE,
                           partition = fx$partition, baseline = fx$baseline)

  expect_s3_class(tab, "data.table")
  expect_identical(nrow(tab), length(fx$models))
  expect_identical(tab$K, sapply(fx$models, function(m) ref_theta_phi(m)$K))
  expect_true(all(c("K", "R2_SE", "R2_chisq", "R2_dev",
                    "R2_SE_macro", "R2_chisq_macro", "R2_dev_macro",
                    "Z_SE", "pval_SE", "Z_chisq", "pval_chisq",
                    "Z_dev", "pval_dev") %in% names(tab)))

  for (i in seq_along(fx$models)) {
    m <- fx$models[[i]]
    se <- optop_index_se(m, fx$dtm, fx$partition, fx$baseline,
                         macro = TRUE, ztest = TRUE)
    x2 <- optop_index_chisq(m, fx$dtm, fx$partition, fx$baseline,
                            macro = TRUE, ztest = TRUE)
    dv <- optop_index_deviance(m, fx$dtm, fx$partition, fx$baseline,
                               macro = TRUE, ztest = TRUE)
    expect_equal(tab$R2_SE[i], se$r2, tolerance = 1e-12)
    expect_equal(tab$R2_chisq[i], x2$r2, tolerance = 1e-12)
    expect_equal(tab$R2_dev[i], dv$r2, tolerance = 1e-12)
    expect_equal(tab$R2_SE_macro[i], se$r2_macro, tolerance = 1e-12)
    expect_equal(tab$R2_chisq_macro[i], x2$r2_macro, tolerance = 1e-12)
    expect_equal(tab$R2_dev_macro[i], dv$r2_macro, tolerance = 1e-12)
    expect_equal(tab$Z_SE[i], se$ztest$z, tolerance = 1e-12)
    expect_equal(tab$pval_SE[i], se$ztest$pval, tolerance = 1e-12)
    expect_equal(tab$Z_chisq[i], x2$ztest$z, tolerance = 1e-12)
    expect_equal(tab$Z_dev[i], dv$ztest$z, tolerance = 1e-12)
  }
})

test_that("optop_index_table reproduces the single-model word indices", {
  fx <- optop_test_fixture()
  tab <- optop_index_table(fx$models, fx$dtm,
                           metrics = c("se", "chisq", "deviance"),
                           level = "word",
                           partition = fx$partition, baseline = fx$baseline)

  expect_true(all(c("K", "R2_SE_micro_word", "R2_SE_macro_word",
                    "R2_chisq_micro_word", "R2_chisq_macro_word",
                    "R2_dev_micro_word", "R2_dev_macro_word") %in% names(tab)))

  for (i in seq_along(fx$models)) {
    m <- fx$models[[i]]
    se <- optop_index_se(m, fx$dtm, fx$partition, fx$baseline, level = "word")
    x2 <- optop_index_chisq(m, fx$dtm, fx$partition, fx$baseline, level = "word")
    dv <- optop_index_deviance(m, fx$dtm, fx$partition, fx$baseline, level = "word")
    expect_equal(tab$R2_SE_micro_word[i], se$r2_micro_word, tolerance = 1e-12)
    expect_equal(tab$R2_SE_macro_word[i], se$r2_macro_word, tolerance = 1e-12)
    expect_equal(tab$R2_chisq_micro_word[i], x2$r2_micro_word, tolerance = 1e-12)
    expect_equal(tab$R2_chisq_macro_word[i], x2$r2_macro_word, tolerance = 1e-12)
    expect_equal(tab$R2_dev_micro_word[i], dv$r2_micro_word, tolerance = 1e-12)
    expect_equal(tab$R2_dev_macro_word[i], dv$r2_macro_word, tolerance = 1e-12)
  }
})

test_that("optop_index_table returns its result visibly", {
  fx <- optop_test_fixture()
  vis <- withVisible(optop_index_table(fx$models[1], fx$dtm, metrics = "se",
                                       partition = fx$partition,
                                       baseline = fx$baseline))
  expect_true(vis$visible)
  expect_s3_class(vis$value, "data.table")
})

test_that("metric subsets control which columns appear", {
  fx <- optop_test_fixture()
  tab <- optop_index_table(fx$models[1:2], fx$dtm, metrics = "se",
                           partition = fx$partition, baseline = fx$baseline)
  expect_named(tab, c("K", "R2_SE"))

  tab2 <- optop_index_table(fx$models[1:2], fx$dtm,
                            metrics = c("chisq", "deviance"),
                            partition = fx$partition, baseline = fx$baseline)
  expect_named(tab2, c("K", "R2_chisq", "R2_dev"))
})

test_that("internally computed partition/baseline match precomputed ones", {
  fx <- optop_test_fixture()
  tab_pre <- optop_index_table(fx$models, fx$dtm, metrics = "se",
                               partition = fx$partition,
                               baseline = fx$baseline)
  tab_int <- optop_index_table(fx$models, fx$dtm, metrics = "se", c = 5)
  expect_equal(tab_int$R2_SE, tab_pre$R2_SE, tolerance = 1e-12)
})

test_that("reopt = 'se' is routed to the SE metric only", {
  withr::local_options(lifecycle_verbosity = "quiet")
  fx <- optop_test_fixture()
  tab_none <- optop_index_table(fx$models, fx$dtm,
                                metrics = c("se", "chisq"), reopt = "none",
                                partition = fx$partition, baseline = fx$baseline)
  tab_se <- optop_index_table(fx$models, fx$dtm,
                              metrics = c("se", "chisq"), reopt = "se",
                              partition = fx$partition, baseline = fx$baseline)
  ref <- sapply(fx$models, function(m) {
    optop_index_se(m, fx$dtm, fx$partition, fx$baseline, reopt = "se")$r2
  })
  expect_equal(tab_se$R2_SE, ref, tolerance = 1e-12)
  expect_equal(tab_se$R2_chisq, tab_none$R2_chisq, tolerance = 1e-12)
})
