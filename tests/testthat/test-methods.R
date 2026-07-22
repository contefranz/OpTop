# The classed results and their print/plot methods: every return keeps its
# plain list (or data.table) structure so dollar-access is untouched, the
# print headers are stable (snapshots), and the holdout plot builds both
# views without touching a display device.

test_that("results carry their classes without changing structure", {
  fx <- optop_test_fixture()

  expect_s3_class(fx$partition, "optop_partition")
  expect_true(is.list(fx$partition))

  res <- suppressMessages(
    optop_index_deviance(fx$models[[2]], fx$dtm, fx$partition, fx$baseline,
                         macro = TRUE)
  )
  expect_s3_class(res, "optop_index")
  expect_true(is.list(res))
  expect_true(all(c("r2", "r2_macro", "r2_doc", "d_model", "d_null",
                    "K", "metric") %in% names(res)))

  word <- suppressMessages(
    optop_index_deviance(fx$models[[2]], fx$dtm, fx$partition, fx$baseline,
                         level = "word")
  )
  expect_s3_class(word, "optop_index")
  expect_true("r2_word" %in% names(word))

  tab <- suppressMessages(
    optop_index_table(fx$models, fx$dtm, metrics = "deviance",
                      partition = fx$partition, baseline = fx$baseline)
  )
  expect_s3_class(tab, "optop_index_table")
  expect_s3_class(tab, "data.table")
  # data.table semantics survive the subclass
  expect_identical(nrow(tab[K > 0]), nrow(tab))
  expect_false(inherits(as.data.frame(tab), "optop_index_table"))
})

test_that("the print headers are stable", {
  fx <- optop_test_fixture()

  expect_snapshot(print(fx$partition))

  res <- suppressMessages(
    optop_index_deviance(fx$models[[2]], fx$dtm, fx$partition, fx$baseline,
                         macro = TRUE)
  )
  expect_snapshot(print(res))

  word <- suppressMessages(
    optop_index_deviance(fx$models[[2]], fx$dtm, fx$partition, fx$baseline,
                         level = "word")
  )
  expect_snapshot(print(word))

  ret <- print(res)
  expect_identical(class(ret), class(res))
})

test_that("holdout and moment-test print through their headers", {
  hx <- optop_holdout_fixture()
  ho <- suppressMessages(
    optop_index_holdout(hx$models, hx$ev, hx$baseline, c = 1,
                        metrics = "deviance")
  )
  expect_message(print(ho), "held-out discrepancy indices")

  mt <- suppressMessages(
    optop_moment_test(hx$models, hx$tr, hx$ev, type = "strata")
  )
  expect_message(print(mt), "Wald")
})

test_that("plot.optop_holdout builds both views", {
  hx <- optop_holdout_fixture()
  ho <- suppressMessages(
    optop_index_holdout(hx$models, hx$ev, hx$baseline, c = 1,
                        metrics = "deviance")
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  p_macro <- plot(ho)
  expect_s3_class(p_macro, "ggplot")
  p_gain <- plot(ho, which = "gains", epsilon = 0.05)
  expect_s3_class(p_gain, "ggplot")
  expect_error(plot(ho, which = "gains", metric = "chisq"),
               "not evaluated")
})
