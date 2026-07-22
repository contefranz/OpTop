# Generates the data behind vignettes/OpTop.Rmd.
#
# The fitted LDA grid is heavy (~170 MB) and stays in data-raw/ (git-ignored;
# re-running this script regenerates it, but the full-corpus fits are unseeded,
# so their exact numbers can shift). The held-out training grid below is seeded.
# Only the small derived tables — plus the captured
# console output of the optop_select() calls, used by the vignette as a
# fallback when the model cache is absent — ship with the package in
# inst/extdata/optop-vignette-results.rds.
#
# CRAN mechanics (deliberate design): data-raw/ is .Rbuildignore'd, so the
# model cache never enters the tarball. The maintainer's R CMD build runs the
# vignette LIVE against the cache and ships that HTML in inst/doc; CRAN's
# "re-building of vignette outputs" re-knits inside the pruned tarball, where
# have_models is FALSE, so the vignette takes the fast replay path (the
# captured console streams stored in the bundle) and rebuilds in seconds with
# identical numbers. Should CRAN's rebuild environment ever become a problem,
# the fallback plan is the rOpenSci "precompiled vignette" pattern
# (OpTop.Rmd.orig knitted locally into a fully static OpTop.Rmd).
#
# Run from the package root: Rscript data-raw/vignette-data.R

library(OpTop)

ncores <- 6

## Part 1 -- U.S. Presidential Inaugural Address corpus -----------------------

# Tokenize
toks <- quanteda::data_corpus_inaugural |>
  quanteda::tokens(remove_punct = TRUE,
                   remove_symbols = TRUE,
                   remove_numbers = TRUE) |>
  quanteda::tokens_tolower() |>
  quanteda::tokens_remove(quanteda::stopwords())

# Create DFM (raw counts - DO NOT weight, OpTop expects counts not proportions)
mydfm        <- quanteda::dfm(toks)
mydfm_sub    <- quanteda::dfm_trim(mydfm, min_termfreq = 5)
weighted_dfm <- quanteda::dfm_weight(mydfm_sub, scheme = "prop")

K_grid <- seq(2, 200, by = 2)

models_path <- file.path("data-raw", "VEM_models_inaugural.rds")
if (file.exists(models_path)) {
  message("Reusing cached ", models_path)
  VEM_models <- readRDS(models_path)
} else {
  # Estimate VEM models
  VEM_models <- parallel::mclapply(
    K_grid,
    function(k) {
      message("Estimating LDA with K = ", k)
      topicmodels::LDA(x = mydfm_sub, k = k)
    },
    mc.cores = ncores
  )
  saveRDS(VEM_models, models_path)
}

alpha <- 0.01

# capture the cli console stream of a call: the vignette replays it verbatim
# when the model cache is not available to evaluate the call live
capture_cli <- function(expr) {
  msgs <- character()
  out <- withCallingHandlers(
    expr,
    message = function(m) {
      msgs <<- c(msgs, cli::ansi_strip(conditionMessage(m)))
      invokeRestart("muffleMessage")
    }
  )
  msgs <- gsub("\r", "", msgs)
  msgs <- unlist(strsplit(msgs, "\n", fixed = TRUE))
  list(value = out, console = msgs[nzchar(msgs)])
}

seq_run <- capture_cli(
  optop_select(topic_models = VEM_models,
               weighted_dfm = weighted_dfm,
               selection = "sequential",
               q = 0.95,
               alpha = alpha,
               do_plot = FALSE)
)
min_run <- capture_cli(
  optop_select(topic_models = VEM_models,
               weighted_dfm = weighted_dfm,
               selection = "min",
               q = 0.95,
               alpha = alpha,
               do_plot = FALSE)
)

# calibrated runs: bootstrap (captured for the vignette, seeded so a live
# vignette build reproduces it exactly) and the closed-form moment matching
doc_lens <- quanteda::ntoken(mydfm_sub)
boot_seed <- 20260703
cal_time <- system.time(
  cal_run <- capture_cli(
    optop_select(topic_models = VEM_models,
                 weighted_dfm = weighted_dfm,
                 q = 0.95,
                 alpha = alpha,
                 calibrate = "bootstrap",
                 n_boot = 200,
                 doc_lengths = doc_lens,
                 seed = boot_seed,
                 do_plot = FALSE)
  )
)
mm_run <- capture_cli(
  optop_select(topic_models = VEM_models,
               weighted_dfm = weighted_dfm,
               q = 0.95,
               alpha = alpha,
               calibrate = "moment",
               doc_lengths = doc_lens,
               do_plot = FALSE)
)

# the returned table does not depend on the selection rule, so one call per q
# is enough for the sweep; the picks mirror the documented rules
q_sweep <- c(0.80, 0.90, 0.95, 0.99)

pick_sequential <- function(tab, alpha) {
  k <- tab$topic[tab$pval > alpha][1]
  if (is.na(k)) tab$topic[which.min(tab$OpTop)] else k
}
pick_min <- function(tab) tab$topic[which.min(tab$OpTop)]

chi_by_q <- lapply(q_sweep, function(q) {
  if (q == 0.95) return(seq_run$value)
  optop_select(VEM_models, weighted_dfm, q = q, alpha = alpha,
               do_plot = FALSE, verbose = TRUE)
})
names(chi_by_q) <- sprintf("q%.2f", q_sweep)

picks <- data.frame(
  q = q_sweep,
  sequential = vapply(chi_by_q, pick_sequential, numeric(1), alpha = alpha),
  min = vapply(chi_by_q, pick_min, numeric(1))
)

## Part 1b -- discrepancy indices on the inaugural corpus ---------------------

dtm_counts <- methods::as(mydfm_sub, "CsparseMatrix")
partition <- optop_make_partition(VEM_models, dtm_counts, c = 1)
baseline <- optop_make_baseline(dtm_counts)
index_tab <- optop_index_table(VEM_models, dtm_counts,
                               metrics = c("se", "chisq", "deviance"),
                               partition = partition, baseline = baseline,
                               macro = TRUE)

# word-level snapshot at the sequential pick (q = 0.95): best and worst words
k_star <- picks$sequential[picks$q == 0.95]
m_star <- VEM_models[[which(K_grid == k_star)]]
word_idx <- optop_index_chisq(m_star, dtm_counts, partition, baseline,
                              level = "word")
r2w <- sort(word_idx$r2_word)
word_snapshot <- list(k = k_star,
                      worst = utils::head(r2w, 10),
                      best = utils::tail(r2w, 10),
                      micro = word_idx$r2_micro_word,
                      macro = word_idx$r2_macro_word)

# null-discrepancy floor demonstration at the sequential pick: the corpus
# genuinely contains a structurally collapsed speech whose support is swept
# entirely into the min bin; the sparse kernels compute its baseline
# discrepancy as an exact zero (complement sums), so the demo shows the
# reported, principled exclusion on real data
dev_floor <- optop_index_deviance(m_star, dtm_counts, partition, baseline,
                                  macro = TRUE)
excl <- is.na(dev_floor$r2_doc) & is.finite(dev_floor$d_null)
null_floor <- list(
  k = k_star,
  n_excluded = dev_floor$n_null_excluded,
  share = dev_floor$null_excluded_share,
  docs = rownames(dtm_counts)[excl],
  lengths = unname(Matrix::rowSums(dtm_counts)[excl]),
  d_null_excluded = unname(dev_floor$d_null[excl]),
  macro_floor = dev_floor$r2_macro,
  micro_floor = dev_floor$r2
)

# WarpLDA demonstration: the same corpus through a different engine, wrapped
# by optop_warplda(); only compact numbers ship with the bundle
set.seed(20260722)
warp_lda <- text2vec::LDA$new(n_topics = k_star)
warp_theta <- warp_lda$fit_transform(dtm_counts, n_iter = 500,
                                     progressbar = FALSE)
warp_fit <- optop_warplda(warp_lda, warp_theta)
warp_part <- optop_make_partition(list(warp_fit), dtm_counts, c = 1)
warp_dev <- optop_index_deviance(warp_fit, dtm_counts, warp_part, baseline,
                                 macro = TRUE)
warplda_demo <- list(k = k_star,
                     micro = warp_dev$r2,
                     macro = warp_dev$r2_macro,
                     micro_vem = dev_floor$r2,
                     macro_vem = dev_floor$r2_macro)

## Part 1c -- held-out validation on the inaugural corpus ----------------------

# Seeded split: about 25% of the speeches are held out for evaluation and the
# grid is refitted on the training split alone, so the evaluation documents
# play no role in estimation.
split_seed <- 20260709
set.seed(split_seed)
n_docs <- quanteda::ndoc(mydfm_sub)
eval_idx <- sort(sample(n_docs, round(0.25 * n_docs)))

dfm_eval <- mydfm_sub[eval_idx, ]
# drop the features that never occur in the training split; words that occur
# only in held-out speeches are routed to the synthetic out-of-vocabulary bin
# by the holdout layer
dfm_train <- quanteda::dfm_trim(mydfm_sub[-eval_idx, ], min_termfreq = 1)

K_grid_train <- seq(2, 60, by = 2)
train_path <- file.path("data-raw", "VEM_models_inaugural_train.rds")
if (file.exists(train_path)) {
  message("Reusing cached ", train_path)
  VEM_models_train <- readRDS(train_path)
} else {
  VEM_models_train <- parallel::mclapply(
    K_grid_train,
    function(k) {
      message("Estimating training-split LDA with K = ", k)
      topicmodels::LDA(x = dfm_train, k = k, control = list(seed = 500 + k))
    },
    mc.cores = ncores
  )
  saveRDS(VEM_models_train, train_path)
}

dtm_train <- methods::as(dfm_train, "CsparseMatrix")
dtm_eval <- methods::as(dfm_eval, "CsparseMatrix")

train_K_actual <- vapply(VEM_models_train, function(model) {
  OpTop:::optop_as_theta_phi(model)$K
}, numeric(1))
if (!identical(train_K_actual, as.numeric(K_grid_train))) {
  stop("cached training-split models do not match K_grid_train; remove and rebuild the cache")
}
if (!identical(OpTop:::optop_as_theta_phi(VEM_models_train[[1L]])$terms,
               colnames(dtm_train))) {
  stop("cached training-split models do not match the current training vocabulary")
}

# out-of-vocabulary share of the evaluation sample, quoted in the vignette
ev_counts <- Matrix::colSums(dtm_eval)
oov_feats <- names(ev_counts)[ev_counts > 0 &
                                !(names(ev_counts) %in% colnames(dtm_train))]
oov_share <- sum(ev_counts[oov_feats]) / sum(ev_counts)

base_train <- optop_make_baseline(dtm_train)
ho <- optop_index_holdout(VEM_models_train, dtm_eval, base_train)
gt <- optop_gain_table(ho, metric = "deviance", epsilon = 0.01)

k_hat_ho <- gt$k_hat
selection_status <- if (is.na(k_hat_ho)) "no_adequate_step" else "selected"
k_diagnostic <- k_hat_ho
if (is.na(k_diagnostic)) {
  # This fallback is solely a topic count at which to illustrate the moment
  # diagnostic. It is not an epsilon-adequacy selection.
  summ_dev <- ho$summary[ho$summary$metric == "deviance", ]
  k_diagnostic <- summ_dev$K[which.max(summ_dev$macro)]
}
m_star_train <- VEM_models_train[[which(K_grid_train == k_diagnostic)]]
mt <- optop_moment_test(list(m_star_train), dtm_eval, dtm_train,
                        type = "strata", bins = 5)

## Part 1d -- V-fold cross-fitting on the inaugural corpus --------------------

# Coarse grid, V = 5: 25 VEM fits, cached per (fold signature, K) so the
# bundle regenerates without refitting when folds and grid are unchanged.
cf_seed <- 20260722
K_grid_cf <- c(10L, 20L, 30L, 40L, 50L)
cf_cache_path <- file.path("data-raw", "VEM_models_crossfit.rds")
cf_cache <- if (file.exists(cf_cache_path)) {
  readRDS(cf_cache_path)
} else {
  list()
}
fit_cf <- function(dtm_train, k) {
  key <- paste(k, nrow(dtm_train), sum(dtm_train), sep = "|")
  if (is.null(cf_cache[[key]])) {
    message("  fitting crossfit VEM K = ", k, " on ", nrow(dtm_train),
            " documents")
    cf_cache[[key]] <<- topicmodels::LDA(
      quanteda::convert(quanteda::as.dfm(as.matrix(dtm_train)),
                        to = "topicmodels"),
      k = k, method = "VEM", control = list(seed = 3000 + k)
    )
  }
  cf_cache[[key]]
}
cf <- optop_crossfit(dtm_counts, K = K_grid_cf, fit_fun = fit_cf, V = 5,
                     metrics = "deviance", c = 1, seed = cf_seed,
                     verbose = TRUE)
saveRDS(cf_cache, cf_cache_path)
crossfit_demo <- list(
  summary = cf$summary,
  gains = cf$gains$deviance$gains,
  k_hat = cf$gains$deviance$k_hat,
  epsilon = cf$gains$deviance$epsilon,
  V = cf$V,
  K = cf$K,
  seed = cf_seed,
  n_docs = length(cf$folds)
)

## Part 2 -- simulation with known true K -------------------------------------

set.seed(20260703)
true_K <- 10L
J_sim <- 100L
W_sim <- 1500L

rdirichlet <- function(n, alpha) {
  g <- matrix(rgamma(n * length(alpha), shape = alpha, rate = 1),
              nrow = n, byrow = TRUE)
  g / rowSums(g)
}

Theta_true <- rdirichlet(J_sim, rep(0.3, true_K))          # J x K
Beta_true <- rdirichlet(true_K, rep(0.05, W_sim))          # K x W, rows sum to 1
colnames(Beta_true) <- sprintf("w%04d", seq_len(W_sim))
doc_len <- sample(500:1500, J_sim, replace = TRUE)

sim_corpus <- sim_dfm(DTW = Theta_true, TWW = Beta_true,
                      doc_length = doc_len, seed = 42)

sim_grid <- seq(2, 20, by = 2)
sim_models <- parallel::mclapply(
  sim_grid,
  function(k) {
    message("Estimating simulated LDA with K = ", k)
    topicmodels::LDA(x = sim_corpus, k = k, control = list(seed = 2000 + k))
  },
  mc.cores = ncores
)

sim_weighted <- quanteda::dfm_weight(sim_corpus, scheme = "prop")
sim_chi <- optop_select(sim_models, sim_weighted, q = 0.95, alpha = alpha,
                        do_plot = FALSE, verbose = TRUE)
sim_picks <- data.frame(
  sequential = pick_sequential(sim_chi, alpha),
  min = pick_min(sim_chi)
)

sim_dtm <- methods::as(sim_corpus, "CsparseMatrix")
sim_partition <- optop_make_partition(sim_models, sim_dtm, c = 1)
sim_baseline <- optop_make_baseline(sim_dtm)
sim_index_tab <- optop_index_table(sim_models, sim_dtm,
                                   metrics = c("se", "chisq", "deviance"),
                                   partition = sim_partition,
                                   baseline = sim_baseline,
                                   macro = TRUE)

## Bundle ----------------------------------------------------------------------

bundle <- list(
  meta = list(
    generated = Sys.time(),
    K_grid = K_grid,
    q_sweep = q_sweep,
    alpha = alpha,
    ndoc = quanteda::ndoc(mydfm_sub),
    nfeat_full = quanteda::nfeat(mydfm),
    nfeat_trimmed = quanteda::nfeat(mydfm_sub),
    true_K = true_K,
    J_sim = J_sim,
    W_sim = W_sim,
    sim_grid = sim_grid,
    boot_seed = boot_seed,
    cal_seconds = unname(cal_time[3])
  ),
  inaugural = list(
    chi_by_q = chi_by_q,
    picks = picks,
    console_seq = seq_run$console,
    console_min = min_run$console,
    chi_cal = cal_run$value,
    console_cal = cal_run$console,
    chi_cal_mm = mm_run$value,
    console_cal_mm = mm_run$console,
    index_table = index_tab,
    word_snapshot = word_snapshot,
    null_floor = null_floor,
    warplda = warplda_demo,
    crossfit = crossfit_demo,
    holdout = list(
      summary = ho$summary,
      gains = gt$gains,
      k_hat = k_hat_ho,
      selection_status = selection_status,
      k_diagnostic = k_diagnostic,
      epsilon = gt$epsilon,
      alpha_gain = gt$alpha,
      moment_summary = mt$summary,
      moment_marginal = mt$moments[[1L]]$marginal,
      split = list(seed = split_seed,
                   n_train = n_docs - length(eval_idx),
                   n_eval = length(eval_idx),
                   K_grid = K_grid_train,
                   nfeat_train = ncol(dtm_train),
                   n_oov = length(oov_feats),
                   oov_share = oov_share)
    )
  ),
  sim = list(
    chi = sim_chi,
    picks = sim_picks,
    index_table = sim_index_tab
  )
)

out_path <- file.path("inst", "extdata", "optop-vignette-results.rds")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(bundle, out_path, compress = "xz")
message("Saved ", out_path, " (", round(file.size(out_path) / 1024), " KB)")
