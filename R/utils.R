# topicmodels adapters ----

# handy function to check LDA_VEM class
is.LDA_VEM <- function(x) {
  "LDA_VEM" %in% class(x)
}
#' @keywords internal


#' @keywords internal
optop_as_theta_phi <- function(model) UseMethod("optop_as_theta_phi")

#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA
optop_as_theta_phi.LDA <- function(model) {
  p <- topicmodels::posterior(model)
  list(theta = p$topics,  # J x K
       phi   = p$terms,   # K x W
       K     = ncol(p$topics))
}

#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA_VEM
optop_as_theta_phi.LDA_VEM <- function(model) optop_as_theta_phi.LDA(model)

#' @keywords internal
#' @exportS3Method optop_as_theta_phi LDA_Gibbs
optop_as_theta_phi.LDA_Gibbs <- function(model) optop_as_theta_phi.LDA(model)

# Fallback for other topicmodels classes (defensive)
#' @keywords internal
#' @exportS3Method optop_as_theta_phi TopicModel
optop_as_theta_phi.TopicModel <- function(model) {
  if (inherits(model, "LDA")) return(optop_as_theta_phi.LDA(model))
  stop("Unsupported 'topicmodels' class: ", paste(class(model), collapse = "/"))
}

#' @keywords internal
# text2vec::LDA (placeholder — second priority)
optop_as_theta_phi.LDA_t2v <- function(model) {
  # TODO: map model$doc_topic_distr -> theta; model$topic_word_distribution -> phi (rows sum to 1)
  stop("text2vec adapter not implemented yet.")
}

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


# this declare my personal ggplot theme
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
