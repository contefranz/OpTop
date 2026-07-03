#' Harmonized Rare-word Partition for OpTop indices
#'
#' Compute the fixed, document-specific rare-word sets \eqn{C_j^*} used to
#' evaluate OpTop goodness-of-fit indices on a **common support** across a grid
#' of topic models. A word \eqn{w} is marked rare in document \eqn{j} if
#' \eqn{\min_K i^{(K)}_{jw} < \tau_j}, where \eqn{i^{(K)}_{jw}} is the fitted
#' probability of word \eqn{w} under model \eqn{K} and
#' \eqn{\tau_j = c / L_j} with \eqn{L_j} the document length.
#'
#' @param models A list of fitted topic models spanning a grid of \eqn{K}.
#'   Currently supports \code{topicmodels::LDA} objects (VEM or Gibbs).
#'   Additional adapters may be added in future versions.
#' @param dtm A document–term matrix of **counts** with rows = documents and
#'   columns = vocabulary. Recommended class is \code{Matrix::dgCMatrix}.
#' @param c Positive scalar controlling the per-document threshold
#'   \eqn{\tau_j = c / L_j}. Larger values classify more words as rare.
#'   Default: \code{5}.
#' @param block Integer; number of terms to process per block when multiplying
#'   \eqn{\Theta \Phi} to control memory usage for large vocabularies.
#'   Default: \code{5000}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{rare_mask}: logical matrix \code{J x W}; \code{TRUE} where
#'   word \code{w} is rare in document \code{j} (belongs to \eqn{C_j^*}).
#'   \item \code{L}: numeric vector of length \code{J} with document lengths
#'   \eqn{L_j = \sum_w N_{jw}}.
#' }
#'
#' @details
#' OpTop indices (SE, Pearson \eqn{\chi^2}, and Deviance) are computed on the
#' fixed support \eqn{\{w \not\in C_j^*\} \cup \{\min\}} for each document, where
#' all words in \eqn{C_j^*} are collapsed into a single "\code{min}" bin. Using
#' a harmonized partition ensures comparability across \eqn{K} and prevents
#' support-driven artifacts in goodness-of-fit curves.
#'
#' @section Computational notes:
#' The function forms \eqn{\min_K (\Theta^{(K)} \Phi^{(K)})} without storing a
#' full \code{J x W x |K|} tensor by updating a running minimum in vocabulary
#' blocks of size \code{block}. For large corpora, keep \code{dtm} sparse and
#' increase \code{block} cautiously.
#'
#' @seealso
#' \code{\link{optop_make_baseline}},
#' \code{\link{optop_index_se}},
#' \code{\link{optop_index_chisq}},
#' \code{\link{optop_index_deviance}},
#' \code{\link{optop_index_table}}
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' library(topicmodels)
#'
#' # toy sparse DTM (J=50, W=1000)
#' set.seed(1)
#' ij <- cbind(
#'   sample(50, 5000, TRUE),
#'   sample(1000, 5000, TRUE)
#' )
#' x  <- sample(1:3, 5000, TRUE)
#' dtm <- sparseMatrix(i = ij[,1], j = ij[,2], x = x, dims = c(50, 1000))
#'
#' # fit a small grid of LDA models
#' m5  <- LDA(dtm, k = 5,  method = "VEM", control = list(seed = 42))
#' m10 <- LDA(dtm, k = 10, method = "VEM", control = list(seed = 42))
#'
#' part <- optop_make_partition(list(m5, m10), dtm, c = 5)
#' str(part$rare_mask)
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' @export
optop_make_partition <- function(models, dtm, c = 5, block = 5000) {
  stopifnot(length(models) >= 1)
  L <- Matrix::rowSums(dtm)                           # L_j
  tau <- c / pmax(L, 1)                               # τ_j
  # initialize running min of fitted probs i_min(j,w) = min_K i^K_{jw}
  theta_phi <- lapply(models, optop_as_theta_phi)
  J <- nrow(theta_phi[[1]]$theta); W <- ncol(theta_phi[[1]]$phi)
  i_min <- matrix(Inf, nrow = J, ncol = W)
  
  # label i_min so that rare_mask inherits the vocabulary and document names
  colnames(i_min) <- colnames(theta_phi[[1]]$phi)
  rownames(i_min) <- rownames(theta_phi[[1]]$theta)

  for (tp in theta_phi) {
    theta <- tp$theta; phi <- tp$phi
    for (start in seq(1, W, by = block)) {
      end <- min(start + block - 1, W)
      I <- theta %*% phi[, start:end, drop = FALSE]   # J x block
      i_min[, start:end] <- pmin(i_min[, start:end], I)
    }
  }
  # rare if min_K i^K_{jw} < tau_j
  rare_mask <- i_min < tau
  list(rare_mask = rare_mask, L = L)
}

#' Global Corpus Baseline for OpTop Indices
#'
#' Compute the global word distribution \eqn{\pi_{\mathrm{glob}}} from a corpus
#' of counts. This baseline defines the "no-topics" model used to normalize
#' OpTop indices via \eqn{B_{jw} = L_j \, \pi_{\mathrm{glob}}(w)} and
#' \eqn{B_{j,\min} = L_j \sum_{w \in C_j^*} \pi_{\mathrm{glob}}(w)}.
#'
#' @param dtm A document–term matrix of **counts** with rows = documents and
#'   columns = vocabulary. Recommended class is \code{Matrix::dgCMatrix}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{pi_glob}: numeric vector of length \code{W} with
#'   \eqn{\pi_{\mathrm{glob}}(w) = N_{\cdot w} / N_{\cdot \cdot}}, i.e., the
#'   corpus-level term probabilities.
#' }
#'
#' @details
#' The returned \code{pi_glob} is used by the exported index functions to form
#' per-document baseline counts on the same evaluation support defined by
#' \code{optop_make_partition()}. Keeping this baseline fixed across \eqn{K}
#' enables interpretable, regression-style \(R^2\) measures.
#'
#' @seealso
#' \code{\link{optop_make_partition}},
#' \code{\link{optop_index_se}},
#' \code{\link{optop_index_chisq}},
#' \code{\link{optop_index_deviance}},
#' \code{\link{optop_index_table}}
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' set.seed(1)
#' dtm <- rsparsematrix(100, 5000, density = 0.002, rand.x = function(n) rpois(n, 2))
#' base <- optop_make_baseline(dtm)
#' sum(base$pi_glob)  # ~ 1
#' }
#' 
#' @export
optop_make_baseline <- function(dtm) {
  N_tot <- Matrix::colSums(dtm)
  L_tot <- sum(N_tot)
  pi_glob <- as.numeric(N_tot) / L_tot
  names(pi_glob) <- colnames(dtm)   # <-- add names for safe alignment
  list(pi_glob = pi_glob)
}
