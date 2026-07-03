# topicmodels adapters ----

#' Check for a topicmodels VEM fit
#'
#' Input gate used by the exported functions: `TRUE` when `x` is an
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
#' Adapter generic behind the discrepancy indices: every supported topic
#' model class is reduced to the common contract
#' `list(theta, phi, K)`, with \eqn{\theta}{theta} the
#' \eqn{J \times K}{J x K} document-topic matrix (rows summing to 1),
#' \eqn{\phi}{phi} the \eqn{K \times W}{K x W} topic-word matrix (rows
#' summing to 1), and \eqn{K} the number of topics. Adding support for a new
#' topic-model implementation means adding a method for its class.
#'
#' @param model A fitted topic model.
#'
#' @return A list with elements `theta`, `phi` and `K` as described above.
#'
#' @keywords internal
optop_as_theta_phi <- function(model) UseMethod("optop_as_theta_phi")

#' @describeIn optop_as_theta_phi `topicmodels::LDA` fits: theta and phi
#'   come from `topicmodels::posterior()`.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA
optop_as_theta_phi.LDA <- function(model) {
  p <- topicmodels::posterior(model)
  list(theta = p$topics,  # J x K
       phi   = p$terms,   # K x W
       K     = ncol(p$topics))
}

#' @describeIn optop_as_theta_phi VEM fits delegate to the `LDA` method.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA_VEM
optop_as_theta_phi.LDA_VEM <- function(model) optop_as_theta_phi.LDA(model)

#' @describeIn optop_as_theta_phi Gibbs fits delegate to the `LDA` method.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA_Gibbs
optop_as_theta_phi.LDA_Gibbs <- function(model) optop_as_theta_phi.LDA(model)

#' @describeIn optop_as_theta_phi Defensive fallback for other `topicmodels`
#'   classes: LDA subclasses are accepted, everything else fails with an
#'   informative error.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi TopicModel
optop_as_theta_phi.TopicModel <- function(model) {
  if (inherits(model, "LDA")) return(optop_as_theta_phi.LDA(model))
  stop("Unsupported 'topicmodels' class: ", paste(class(model), collapse = "/"))
}

#' @describeIn optop_as_theta_phi Placeholder for `text2vec` (WarpLDA)
#'   support, registered so the method dispatches and fails informatively
#'   until the adapter is implemented.
#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA_t2v
optop_as_theta_phi.LDA_t2v <- function(model) {
  # TODO: map model$doc_topic_distr -> theta; model$topic_word_distribution -> phi (rows sum to 1)
  stop("text2vec adapter not implemented yet.")
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
  vocabs <- lapply(models, function(m) colnames(topicmodels::posterior(m)$terms))
  vocab  <- Reduce(intersect, vocabs)
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
