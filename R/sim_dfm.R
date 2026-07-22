#' Simulate a document-feature matrix from an LDA specification
#'
#' Draw a synthetic corpus of word counts from the standard LDA generative
#' model: each document \eqn{j} mixes the topic-word distributions with its
#' topic weights, \eqn{p_j = \theta_j^\top \Phi}{p_j = theta_j' Phi}, and its
#' \eqn{L_j} tokens are one multinomial draw from \eqn{p_j}. The number of
#' topics, documents, and vocabulary size are inferred from the inputs. The
#' generator is self-contained (no external simulation dependency) and is the
#' quickest way to build small corpora with a known truth, for examples and
#' for exercises like the vignette's known-truth study.
#'
#' @param DTW A matrix or data.frame of document-topic weights
#'   (\eqn{J \times k}{J x k}), rows summing to 1. Row names, when present,
#'   become the document names of the result. Ignored (except for its
#'   dimensions) when `alpha` is supplied.
#' @param TWW A matrix or data.frame of topic-word weights
#'   (\eqn{k \times W}{k x W}), rows summing to 1. Column names, when
#'   present, become the feature names of the result.
#' @param doc_length Numeric vector with the desired length (total token
#'   count) of each document; a single value is recycled to all \eqn{J}
#'   documents.
#' @param alpha Optional single positive number or vector of \eqn{k}
#'   positive numbers; when supplied, the document-topic weights are drawn
#'   fresh from a \eqn{\mathrm{Dirichlet}(\alpha)}{Dirichlet(alpha)} instead
#'   of using `DTW` (default `NULL`, use `DTW` as given).
#' @param seed Optional integer passed to [set.seed()] for reproducible
#'   draws (default `NULL`).
#'
#' @return A [quanteda::dfm] of simulated word counts with one row per
#'   document.
#'
#' @section Change in 0.16.0:
#' The generator is now implemented inside OpTop rather than delegated to
#' the LDATS package. The simulated distribution is unchanged (LDA
#' Dirichlet-multinomial sampling), but the random draws for a given `seed`
#' differ from those of earlier versions.
#'
#' @examples
#' # a small corpus with known truth: 20 documents, 3 topics, 50 words
#' rdirich <- function(n, k) {
#'   g <- matrix(stats::rgamma(n * k, shape = 1), n, k)
#'   g / rowSums(g)
#' }
#' theta <- rdirich(20, 3)
#' phi <- rdirich(3, 50)
#' colnames(phi) <- sprintf("w%02d", 1:50)
#'
#' sim <- sim_dfm(theta, phi, doc_length = 200, seed = 42)
#' sim
#' quanteda::ntoken(sim)[1:5]
#'
#' @seealso [optimal_topic()], [optop_index_table()]
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' @export
sim_dfm <- function(DTW, TWW, doc_length, alpha = NULL, seed = NULL) {
  if (!is.matrix(DTW) && !is.data.frame(DTW)) {
    stop("DTW must be either a matrix or a data.frame")
  }
  if (!is.matrix(TWW) && !is.data.frame(TWW)) {
    stop("TWW must be either a matrix or a data.frame")
  }
  if (!is.numeric(doc_length) || length(doc_length) == 0L ||
      any(!is.finite(doc_length)) || any(doc_length <= 0)) {
    stop("doc_length must be a vector of positive document lengths")
  }
  theta <- as.matrix(DTW)
  phi <- as.matrix(TWW)
  J <- nrow(theta)
  k <- ncol(theta)
  if (nrow(phi) != k) {
    stop("TWW must have one row per topic: nrow(TWW) == ncol(DTW)")
  }
  if (length(doc_length) == 1L) {
    doc_length <- rep(doc_length, J)
  }
  if (length(doc_length) != J) {
    stop("doc_length must have one entry per document of DTW")
  }
  sizes <- as.integer(round(doc_length))
  doc_names <- rownames(theta)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(alpha)) {
    if (!is.numeric(alpha) || !(length(alpha) %in% c(1L, k)) ||
        any(!is.finite(alpha)) || any(alpha <= 0)) {
      stop("alpha must be a positive number or a vector of k positive numbers")
    }
    # theta_j ~ Dirichlet(alpha) via normalized gamma draws
    shape <- if (length(alpha) == 1L) rep(alpha, k) else alpha
    g <- matrix(stats::rgamma(J * k, shape = rep(shape, each = J)), J, k)
    theta <- g / rowSums(g)
  } else if (any(theta < 0) || any(abs(rowSums(theta) - 1) > 1e-6)) {
    stop("DTW rows must be topic distributions summing to 1")
  }
  if (any(phi < 0) || any(abs(rowSums(phi) - 1) > 1e-6)) {
    stop("TWW rows must be word distributions summing to 1")
  }

  # p_j = theta_j' Phi; one multinomial draw per document
  P <- theta %*% phi
  W <- ncol(phi)
  counts <- matrix(0L, J, W)
  for (j in seq_len(J)) {
    counts[j, ] <- stats::rmultinom(1L, size = sizes[j], prob = P[j, ])
  }
  dimnames(counts) <- list(
    if (is.null(doc_names)) sprintf("text%d", seq_len(J)) else doc_names,
    if (is.null(colnames(phi))) sprintf("feat%d", seq_len(W))
    else colnames(phi)
  )
  quanteda::as.dfm(counts)
}
