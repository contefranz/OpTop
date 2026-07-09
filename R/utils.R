# model adapters ----

#' Check for a topicmodels VEM fit
#'
#' Input gate used by the deprecated functions and by
#' `optimal_topic(selection = "legacy")`: `TRUE` when `x` is an
#' `LDA_VEM` object as returned by `topicmodels::LDA(method = "VEM")`.
#'
#' @param x Any object.
#'
#' @return A logical scalar.
#'
#' @keywords internal
is.LDA_VEM <- function(x) {
  "LDA_VEM" %in% class(x)
}

#' Extract document-topic and topic-word weights from a fitted model
#'
#' Adapter generic behind [optimal_topic()] and the discrepancy indices:
#' every supported topic-model class is reduced to the common contract
#' `list(theta, phi, K, docs, terms)`, with \eqn{\theta}{theta} the
#' \eqn{J \times K}{J x K} document-topic matrix (rows summing to 1),
#' \eqn{\phi}{phi} the \eqn{K \times W}{K x W} topic-word matrix (rows
#' summing to 1), \eqn{K} the number of topics, `docs` the \eqn{J} document
#' identifiers matching the rows of \eqn{\theta}{theta}, and `terms` the
#' \eqn{W} vocabulary entries matching the columns of \eqn{\phi}{phi}.
#' Document identifiers are mandatory: alignment to the data is always by
#' identifier, never positional. Adding support for a new topic-model
#' implementation means adding a method for its class; the contract is
#' enforced by `.optop_validate_theta_phi()`.
#'
#' @param model A fitted topic model.
#'
#' @return A list with elements `theta`, `phi`, `K`, `docs` and `terms` as
#'   described above.
#'
#' @keywords internal
optop_as_theta_phi <- function(model) UseMethod("optop_as_theta_phi")

#' @describeIn optop_as_theta_phi Workhorse for all `topicmodels` fits (LDA
#'   and CTM alike): theta and phi come from `topicmodels::posterior()`,
#'   identifiers from the `documents` and `terms` slots.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi TopicModel
optop_as_theta_phi.TopicModel <- function(model) {
  p <- topicmodels::posterior(model)
  list(theta = p$topics,  # J x K
       phi   = p$terms,   # K x W
       K     = ncol(p$topics),
       docs  = as.character(model@documents),
       terms = as.character(model@terms))
}

#' @describeIn optop_as_theta_phi `topicmodels::LDA` fits delegate to the
#'   `TopicModel` method.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA
optop_as_theta_phi.LDA <- function(model) optop_as_theta_phi.TopicModel(model)

#' @describeIn optop_as_theta_phi VEM fits delegate to the `TopicModel`
#'   method.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA_VEM
optop_as_theta_phi.LDA_VEM <- function(model) optop_as_theta_phi.TopicModel(model)

#' @describeIn optop_as_theta_phi Gibbs fits delegate to the `TopicModel`
#'   method.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA_Gibbs
optop_as_theta_phi.LDA_Gibbs <- function(model) optop_as_theta_phi.TopicModel(model)

#' @describeIn optop_as_theta_phi `topicmodels::CTM` fits delegate to the
#'   `TopicModel` method.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi CTM
optop_as_theta_phi.CTM <- function(model) optop_as_theta_phi.TopicModel(model)

#' @describeIn optop_as_theta_phi CTM VEM fits delegate to the `TopicModel`
#'   method.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi CTM_VEM
optop_as_theta_phi.CTM_VEM <- function(model) optop_as_theta_phi.TopicModel(model)

#' @describeIn optop_as_theta_phi seededlda fits: one method covers
#'   `seededlda::textmodel_lda()`, `seededlda::textmodel_seededlda()` and
#'   `seededlda::textmodel_seqlda()`, since all three constructors return
#'   objects of class `textmodel_lda` exposing `theta` and `phi`.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi textmodel_lda
optop_as_theta_phi.textmodel_lda <- function(model) {
  list(theta = model$theta,
       phi   = model$phi,
       K     = ncol(model$theta),
       docs  = rownames(model$theta),
       terms = colnames(model$phi))
}

#' @describeIn optop_as_theta_phi NLPstudio fits: theta and phi come from the
#'   stored `dtw` (document-topic) and `tww` (topic-word) weight matrices,
#'   identifiers from `doc_ids` and `vocab`. When the weights are not stored,
#'   the adapter recurses into the raw `model_object` if present.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi nlp_topic_fit
optop_as_theta_phi.nlp_topic_fit <- function(model) {
  theta <- model$dtw
  phi <- model$tww
  if (is.null(theta) || is.null(phi)) {
    if (!is.null(model$model_object)) {
      return(optop_as_theta_phi(model$model_object))
    }
    stop(paste("this nlp_topic_fit stores neither its dtw/tww weight",
               "matrices nor the raw model_object, so theta and phi cannot",
               "be recovered"))
  }
  theta <- as.matrix(theta)
  phi <- as.matrix(phi)
  docs <- model$doc_ids
  if (is.null(docs)) {
    docs <- rownames(theta)
  }
  terms <- model$vocab
  if (is.null(terms)) {
    terms <- colnames(phi)
  }
  list(theta = theta,
       phi   = phi,
       K     = ncol(theta),
       docs  = as.character(docs),
       terms = as.character(terms))
}

#' @describeIn optop_as_theta_phi text2vec WarpLDA fits, wrapped by
#'   [optop_warplda()]: theta is the matrix the user kept from
#'   `fit_transform()`, phi is `model$topic_word_distribution`.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi optop_warplda
optop_as_theta_phi.optop_warplda <- function(model) {
  theta <- model$theta
  phi <- model$model$topic_word_distribution
  if (ncol(theta) != nrow(phi)) {
    stop("doc_topic_distr and the WarpLDA fit disagree on the number of topics")
  }
  list(theta = theta,
       phi   = phi,
       K     = ncol(theta),
       docs  = rownames(theta),
       terms = colnames(phi))
}

#' @describeIn optop_as_theta_phi Guard for raw text2vec WarpLDA R6 objects:
#'   they do not retain the document-topic matrix, so they must be wrapped
#'   with [optop_warplda()]. Registered on the R6 class so that raw objects
#'   never fall through to the `topicmodels` `LDA` method (the R6 class
#'   chain contains "LDA" and "TopicModel").
#' @keywords internal
#' @exportS3Method optop_as_theta_phi WarpLDA
optop_as_theta_phi.WarpLDA <- function(model) {
  stop(paste("a raw text2vec WarpLDA object does not retain the",
             "document-topic matrix; wrap it with",
             "optop_warplda(model, doc_topic_distr), where doc_topic_distr",
             "is the matrix returned by model$fit_transform()"))
}

#' @describeIn optop_as_theta_phi Unsupported classes fail with the list of
#'   supported engines.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi default
optop_as_theta_phi.default <- function(model) {
  stop(paste0("unsupported topic model class <",
              paste(class(model), collapse = "/"),
              ">. Supported: topicmodels (LDA VEM/Gibbs, CTM), seededlda ",
              "(textmodel_lda/textmodel_seededlda/textmodel_seqlda), ",
              "NLPstudio (nlp_topic_fit) and text2vec WarpLDA via ",
              "optop_warplda()."))
}

#' @describeIn optop_as_theta_phi Internal container built by the held-out
#'   tools: the adapter contract stored as-is.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi optop_theta_phi
optop_as_theta_phi.optop_theta_phi <- function(model) {
  model[c("theta", "phi", "K", "docs", "terms")]
}

# Internal light container carrying an explicit (theta, phi) pair through
# the adapter interface, used by the held-out construction where theta is
# folded in rather than stored in a fitted object.
.optop_tp <- function(theta, phi) {
  structure(list(theta = theta,
                 phi   = phi,
                 K     = ncol(theta),
                 docs  = rownames(theta),
                 terms = colnames(phi)),
            class = c("optop_theta_phi", "list"))
}

#' Wrap a text2vec WarpLDA fit for OpTop
#'
#' text2vec's WarpLDA sampler returns the document-topic matrix from
#' `fit_transform()` without retaining it in the model object, while the
#' topic-word distribution stays available as an active binding. This helper
#' bundles the two into a light object accepted by every OpTop tool that
#' takes fitted topic models ([optimal_topic()], [optop_make_partition()],
#' the index family, [optop_index_holdout()] and [optop_moment_test()]).
#'
#' @param model A fitted `text2vec::LDA` (WarpLDA) R6 object.
#' @param doc_topic_distr The document-topic probability matrix returned by
#'   `model$fit_transform()` on the corpus under evaluation (rows are
#'   documents and carry the document identifiers, columns are topics, rows
#'   sum to 1).
#'
#' @return An object of class `optop_warplda` holding the pair; its
#'   `optop_as_theta_phi()` method exposes theta and phi with the usual
#'   contract.
#'
#' @examples
#' \dontrun{
#' lda <- text2vec::LDA$new(n_topics = 10)
#' theta <- lda$fit_transform(dtm, n_iter = 1000)
#' fit <- optop_warplda(lda, theta)
#' part <- optop_make_partition(list(fit), dtm)
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' @export
optop_warplda <- function(model, doc_topic_distr) {
  if (!inherits(model, "WarpLDA")) {
    stop("model must be a fitted text2vec WarpLDA object")
  }
  if (!is.matrix(doc_topic_distr) || is.null(rownames(doc_topic_distr))) {
    stop(paste("doc_topic_distr must be the document-topic matrix returned",
               "by model$fit_transform(), with document identifiers as",
               "row names"))
  }
  structure(list(model = model, theta = doc_topic_distr),
            class = c("optop_warplda", "list"))
}

#' Fold document-topic weights in for new documents
#'
#' Adapter generic behind the held-out tools: given a topic model fitted on
#' the training sample and the counts of new documents on the same
#' vocabulary, return the folded-in document-topic weights
#' \eqn{\hat\theta_j}{theta_j} (rows summing to 1, aligned with the rows of
#' `newdata`), holding the trained topic-word distributions fixed. Every
#' method infers the weights from the new document and the trained global
#' objects only, the held-out-document reconstruction target of the paper.
#'
#' @param model A fitted topic model supported by [optop_as_theta_phi()].
#' @param newdata A counts document-term matrix on the training vocabulary.
#' @param ... Passed to the engine's fold-in routine.
#'
#' @return A numeric matrix of document-topic weights, one row per document
#'   of `newdata`.
#'
#' @keywords internal
optop_fold_in <- function(model, newdata, ...) UseMethod("optop_fold_in")

#' @describeIn optop_fold_in topicmodels fits: the fold-in E step of
#'   `topicmodels::posterior()` on the new documents.
#' @keywords internal
#' @exportS3Method optop_fold_in TopicModel
optop_fold_in.TopicModel <- function(model, newdata, ...) {
  x <- methods::as(newdata, "CsparseMatrix")
  trip <- Matrix::summary(x)
  stm <- slam::simple_triplet_matrix(i = trip$i, j = trip$j, v = trip$x,
                                     nrow = nrow(x), ncol = ncol(x),
                                     dimnames = dimnames(x))
  topicmodels::posterior(model, newdata = stm, ...)$topics
}

#' @describeIn optop_fold_in seededlda fits: refit the new documents with
#'   the trained model held fixed (`textmodel_lda(newdata, model = .)`),
#'   the package's supported prediction path.
#' @keywords internal
#' @exportS3Method optop_fold_in textmodel_lda
optop_fold_in.textmodel_lda <- function(model, newdata, ...) {
  refit <- suppressMessages(
    seededlda::textmodel_lda(quanteda::as.dfm(newdata), model = model,
                             verbose = FALSE, ...)
  )
  refit$theta
}

#' @describeIn optop_fold_in WarpLDA wrappers: `model$transform()` on the
#'   new documents.
#' @keywords internal
#' @exportS3Method optop_fold_in optop_warplda
optop_fold_in.optop_warplda <- function(model, newdata, ...) {
  theta <- model$model$transform(newdata, ...)
  rownames(theta) <- rownames(newdata)
  theta
}

#' @describeIn optop_fold_in NLPstudio fits: recurse into the raw backend
#'   fit stored in `model_object`.
#' @keywords internal
#' @exportS3Method optop_fold_in nlp_topic_fit
optop_fold_in.nlp_topic_fit <- function(model, newdata, ...) {
  if (is.null(model$model_object)) {
    stop(paste("this nlp_topic_fit stores no backend model object;",
               "fold-in needs the raw engine fit"))
  }
  optop_fold_in(model$model_object, newdata, ...)
}

#' @describeIn optop_fold_in Unsupported classes fail with the list of
#'   engines that support fold-in.
#' @keywords internal
#' @exportS3Method optop_fold_in default
optop_fold_in.default <- function(model, newdata, ...) {
  stop(paste0("fold-in is not available for class <",
              paste(class(model), collapse = "/"),
              ">. Supported: topicmodels fits, seededlda fits, ",
              "NLPstudio nlp_topic_fit with a stored backend model, and ",
              "text2vec WarpLDA via optop_warplda()."))
}

#' Validate an adapter result
#'
#' Contract enforcement behind [optimal_topic()]: checks that an
#' `optop_as_theta_phi()` result is a coherent
#' `list(theta, phi, K, docs, terms)` before any statistic is computed, so a
#' broken or partial adapter fails with a clear message instead of
#' propagating misaligned numbers. Document identifiers and terms are
#' mandatory and must be unique, and both weight matrices must have rows
#' summing to 1.
#'
#' @param tp A list as returned by `optop_as_theta_phi()`.
#' @param model_class Character; the class vector of the adapted model, used
#'   in the error messages.
#'
#' @return Invisibly `TRUE`; called for its side effect of stopping on a
#'   contract violation.
#'
#' @keywords internal
.optop_validate_theta_phi <- function(tp, model_class) {
  lbl <- paste(model_class, collapse = "/")
  theta <- tp$theta
  phi <- tp$phi
  if (!is.matrix(theta) || !is.numeric(theta) ||
      !is.matrix(phi) || !is.numeric(phi)) {
    stop("adapter for class <", lbl,
         "> must return numeric matrices theta and phi")
  }
  if (ncol(theta) != nrow(phi)) {
    stop("adapter for class <", lbl,
         ">: ncol(theta) and nrow(phi) disagree on the number of topics")
  }
  if (length(tp$K) != 1L || !is.finite(tp$K) ||
      as.integer(tp$K) != ncol(theta)) {
    stop("adapter for class <", lbl,
         ">: K does not match the dimensions of theta and phi")
  }
  docs <- tp$docs
  if (is.null(docs) || length(docs) != nrow(theta) || anyNA(docs)) {
    stop("adapter for class <", lbl, "> must expose one document identifier",
         " per row of theta; positional alignment is never assumed")
  }
  if (anyDuplicated(docs)) {
    stop("model of class <", lbl, "> has duplicated document identifiers;",
         " alignment would be ambiguous")
  }
  terms <- tp$terms
  if (is.null(terms) || length(terms) != ncol(phi) || anyNA(terms)) {
    stop("adapter for class <", lbl,
         "> must expose one term per column of phi")
  }
  if (anyDuplicated(terms)) {
    stop("model of class <", lbl, "> has duplicated terms; feature",
         " alignment would be ambiguous")
  }
  if (any(!is.finite(theta)) || any(!is.finite(phi)) ||
      any(theta < 0) || any(phi < 0)) {
    stop("model of class <", lbl,
         "> has non-finite or negative fitted probabilities")
  }
  if (any(abs(rowSums(theta) - 1) > 1e-6)) {
    stop("model of class <", lbl, ">: the rows of theta do not sum to 1",
         " (document-topic weights must be proportions)")
  }
  if (any(abs(rowSums(phi) - 1) > 1e-6)) {
    stop("model of class <", lbl, ">: the rows of phi do not sum to 1",
         " (topic-word weights must be proportions)")
  }
  invisible(TRUE)
}

#' Align a document-term matrix to a grid of fitted models
#'
#' Remediation helper referenced by the alignment errors of the index
#' functions: subsets and reorders the columns of `dtm` to the intersection
#' of the models' vocabularies, so that counts, partition and baseline can be
#' recomputed on a support common to the whole grid. Accepts a
#' `quanteda::dfm` (converted to `dgCMatrix`) or any matrix with column
#' names.
#'
#' @param dtm A document-term matrix of counts (rows are documents).
#' @param models A list of fitted topic models supported by
#'   `optop_as_theta_phi()`.
#'
#' @return The realigned document-term matrix, with columns ordered as the
#'   common vocabulary.
#'
#' @keywords internal
optop_align_dtm_to_models <- function(dtm, models) {
  # collect model vocabularies
  vocabs <- lapply(models, function(m) optop_as_theta_phi(m)$terms)
  vocab <- Reduce(intersect, vocabs)
  if (length(vocab) == 0L) stop("Empty vocabulary intersection across models.")

  # subset & reorder dtm to model vocab order
  if (inherits(dtm, "dfm")) {
    # quanteda path
    dtm2 <- quanteda::dfm_select(dtm, pattern = vocab, selection = "keep")
    dtm2 <- quanteda::dfm_match(dtm2, features = vocab)  # reorder to vocab
    dtm2 <- methods::as(dtm2, "dgCMatrix")
    colnames(dtm2) <- vocab
  } else {
    idx <- match(vocab, colnames(dtm))
    if (anyNA(idx)) stop("DTM is missing terms present in model vocab.")
    dtm2 <- dtm[, idx, drop = FALSE]
  }
  dtm2
}


# personal ggplot theme shared by the (deprecated) stability plots; not a
# function, so documented here rather than via roxygen
font_size <- 10
theme_OpTop <- ggplot2::theme(
  title = ggplot2::element_text(face = "bold", size = 8),
  axis.title.x = ggplot2::element_text(face = "bold", size = font_size),
  axis.title.y = ggplot2::element_text(face = "bold", size = font_size),
  axis.text.x = ggplot2::element_text(size = font_size),
  axis.text.y = ggplot2::element_text(size = font_size),
  legend.text = ggplot2::element_text(size = 8),
  legend.position = "bottom"
)
