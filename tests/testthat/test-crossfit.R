# V-fold cross-fitting: fold assignment, leakage-free orchestration,
# pooling math against a hand-rolled reference, the guards, and the
# methods. Engine fits are kept tiny (VEM on toy halves); the small-fold
# warning is expected throughout and asserted once.

optop_crossfit_fixture <- function() {
  if (is.null(.optop_test_cache$crossfit)) {
    fx <- optop_test_fixture()
    fits <- new.env(parent = emptyenv())
    fit_rec <- new.env(parent = emptyenv())
    fit_rec$train_rows <- list()
    fit_fun <- function(dtm_train, k) {
      key <- paste(k, nrow(dtm_train), sum(dtm_train),
                   paste(range(Matrix::rowSums(dtm_train)), collapse = "-"),
                   sep = "|")
      fit_rec$train_rows[[length(fit_rec$train_rows) + 1L]] <-
        list(k = k, rows = rownames(dtm_train))
      if (is.null(fits[[key]])) {
        fits[[key]] <- topicmodels::LDA(as.matrix(dtm_train), k = k,
                                        method = "VEM",
                                        control = list(seed = 700 + k))
      }
      fits[[key]]
    }
    cf <- suppressWarnings(suppressMessages(
      optop_crossfit(fx$dtm, K = 2:3, fit_fun = fit_fun, V = 2,
                     metrics = "deviance", c = 1, seed = 99)
    ))
    .optop_test_cache$crossfit <- list(fx = fx, cf = cf, fit_fun = fit_fun,
                                       rec = fit_rec)
  }
  .optop_test_cache$crossfit
}

test_that("fold assignment covers every document once and reproduces under a seed", {
  L <- rexp(101, 1 / 200) + 20
  set.seed(5)
  f1 <- OpTop:::.optop_crossfit_folds(L, 5L, stratify = TRUE)
  set.seed(5)
  f2 <- OpTop:::.optop_crossfit_folds(L, 5L, stratify = TRUE)
  expect_identical(f1, f2)
  expect_identical(sort(unique(f1)), 1:5)
  expect_true(max(table(f1)) - min(table(f1)) <= 1)

  # stratification balances lengths: every fold's mean length sits inside
  # a tight band around the overall mean
  mets <- tapply(L, f1, mean)
  expect_true(all(abs(mets - mean(L)) < stats::sd(L) / 2))

  set.seed(6)
  f3 <- OpTop:::.optop_crossfit_folds(L, 5L, stratify = FALSE)
  expect_identical(sort(unique(f3)), 1:5)
  expect_true(max(table(f3)) - min(table(f3)) <= 1)
})

test_that("cross-fitting never leaks an evaluation document into training", {
  cfx <- optop_crossfit_fixture()
  cf <- cfx$cf
  fx <- cfx$fx

  # every document scored exactly once, out of fold
  expect_identical(sort(names(cf$folds)), sort(rownames(fx$counts)))
  expect_identical(sort(unique(unname(cf$folds))), 1:2)

  calls <- cfx$rec$train_rows
  expect_identical(length(calls), 4L)  # V = 2 folds x |K| = 2
  for (cl in calls) {
    trained_on <- cl$rows
    held_out <- names(cf$folds)[!names(cf$folds) %in% trained_on]
    # the held-out complement of every call is exactly one fold
    expect_identical(sort(held_out),
                     sort(names(cf$folds)[cf$folds ==
                                            cf$folds[held_out[1]]]))
  }

  # scores exist for every document (the toy corpus has no floored docs
  # at c = 1 with these lengths, except any reported exclusions)
  n_scored <- sum(!is.na(cf$scores$deviance[, 1])) +
    cf$summary$n_null_excluded[1]
  expect_identical(as.integer(n_scored), nrow(fx$counts))
})

test_that("the pooled estimates match a hand-rolled two-fold reference", {
  cfx <- optop_crossfit_fixture()
  cf <- cfx$cf
  fx <- cfx$fx

  # rebuild both folds by hand with the same assignment and fitter
  hos <- lapply(1:2, function(v) {
    idx_ev <- which(cf$folds == v)
    dtm_tr <- fx$dtm[-idx_ev, , drop = FALSE]
    dtm_ev <- fx$dtm[idx_ev, , drop = FALSE]
    models <- lapply(2:3, function(k) cfx$fit_fun(dtm_tr, k))
    suppressMessages(
      optop_index_holdout(models, dtm_ev, optop_make_baseline(dtm_tr),
                          c = 1, metrics = "deviance")
    )
  })

  # per-document scores agree fold-wise
  for (v in 1:2) {
    rows_v <- rownames(hos[[v]]$scores$deviance)
    expect_identical(cf$scores$deviance[rows_v, ],
                     hos[[v]]$scores$deviance)
  }

  # pooled macro is the mean over the stacked J+, micro the pooled ratio
  r_all <- do.call(rbind, lapply(hos, function(h) h$scores$deviance))
  d_k <- do.call(rbind, lapply(hos, function(h) h$d_model$deviance))
  d_0 <- do.call(rbind, lapply(hos, function(h) h$d_null$deviance))
  for (i in 1:2) {
    srow <- cf$summary[cf$summary$K == c(2, 3)[i], ]
    expect_equal(srow$macro, mean(r_all[, i], na.rm = TRUE),
                 tolerance = 1e-12)
    valid <- !is.na(r_all[, i])
    expect_equal(srow$micro,
                 1 - sum(d_k[valid, i]) / sum(d_0[valid, i]),
                 tolerance = 1e-12)
    # fold-replicate standard errors from the two fold macros
    m_v <- vapply(hos, function(h) {
      h$summary$macro[h$summary$K == c(2, 3)[i]]
    }, numeric(1))
    expect_equal(srow$macro_se, stats::sd(m_v) / sqrt(2),
                 tolerance = 1e-12)
  }

  # cross-fitted gains: pooled mean delta with the fold-replicate bound
  g <- cf$gains$deviance$gains
  delta <- r_all[, 2] - r_all[, 1]
  expect_equal(g$gain, mean(delta, na.rm = TRUE), tolerance = 1e-12)
  gv <- vapply(1:2, function(v) {
    rows_v <- rownames(hos[[v]]$scores$deviance)
    mean((cf$scores$deviance[rows_v, 2] -
            cf$scores$deviance[rows_v, 1]), na.rm = TRUE)
  }, numeric(1))
  expect_equal(g$se_fold, stats::sd(gv) / sqrt(2), tolerance = 1e-12)
})

test_that("optop_crossfit validates its inputs", {
  fx <- optop_test_fixture()
  dummy_fit <- function(dtm_train, k) stop("never called")

  expect_error(
    optop_crossfit(optop_corpus(fx$dtm), K = 2:3, fit_fun = dummy_fit),
    "optop_corpus"
  )
  expect_error(
    optop_crossfit(fx$dtm, K = c(3, 2), fit_fun = dummy_fit),
    "strictly increasing"
  )
  expect_error(
    optop_crossfit(fx$dtm, K = 2:3, fit_fun = function(x) x),
    "fit_fun"
  )
  expect_error(
    optop_crossfit(fx$dtm, K = 2:3, fit_fun = dummy_fit, V = 40),
    "between 2 and"
  )
  anon <- fx$dtm
  rownames(anon) <- NULL
  expect_error(
    optop_crossfit(anon, K = 2:3, fit_fun = dummy_fit),
    "document names"
  )
  # the small-fold warning fires on toy corpora
  expect_warning(
    suppressMessages(
      optop_crossfit(fx$dtm, K = 2:3, V = 2, seed = 1,
                     metrics = "deviance",
                     fit_fun = function(dtm_train, k) {
                       topicmodels::LDA(as.matrix(dtm_train), k = k,
                                        method = "VEM",
                                        control = list(seed = 1))
                     })
    ),
    "fewer than 30 documents per fold"
  )
})

test_that("a fit_fun returning the wrong K is caught", {
  fx <- optop_test_fixture()
  wrong <- function(dtm_train, k) {
    topicmodels::LDA(as.matrix(dtm_train), k = k + 1, method = "VEM",
                     control = list(seed = 1))
  }
  expect_error(
    suppressWarnings(suppressMessages(
      optop_crossfit(fx$dtm, K = 2:3, fit_fun = wrong, V = 2, seed = 1,
                     metrics = "deviance")
    )),
    "topic counts"
  )
})

test_that("the crossfit result prints and plots", {
  cfx <- optop_crossfit_fixture()
  cf <- cfx$cf

  expect_s3_class(cf, "optop_crossfit")
  expect_message(ret <- print(cf), "cross-fitted")
  expect_identical(class(ret), class(cf))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(plot(cf), "ggplot")
  expect_s3_class(plot(cf, which = "gains", epsilon = 0.5), "ggplot")
  expect_error(plot(cf, which = "gains", metric = "chisq"),
               "not evaluated")
})