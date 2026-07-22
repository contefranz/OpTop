#' Sharded corpus for very large document sets
#'
#' Wrap a document-term matrix, or a collection of document shards of one,
#' into the container the OpTop entry points accept wherever a dtm or a
#' weighted dfm is accepted today. Sharding is how OpTop passes the hard
#' container limit of R sparse matrices: a single `dgCMatrix` (and hence a
#' single quanteda dfm) holds at most \eqn{2^{31} - 1}{2^31 - 1} nonzero
#' entries, a bound no compiler flag lifts. Every OpTop statistic is a sum
#' of independent per-document terms, so evaluating shard by shard and
#' combining is exact, not an approximation: per-document scores are
#' bit-identical to an unsharded run and concatenate; statistics and
#' degrees of freedom add, agreeing with the unsharded totals up to
#' floating-point summation order; and the calibration bootstrap keys
#' every document's random stream by its global index, so each document
#' draws identical replicates whatever the sharding.
#'
#' @param x One of:
#'   * a single dfm or sparse matrix (the trivial one-shard corpus);
#'   * a list of dfms or sparse matrices, the shards, in document order;
#'   * a character vector of file paths, one shard per file, loaded on
#'     demand through `reader` (shards are then held in memory one at a
#'     time).
#' @param reader A function of one path returning a dfm or sparse matrix;
#'   required when `x` is a character vector (for example
#'   `function(p) qs2::qs_read(p)` or [readRDS()]).
#'
#' @details
#' All shards must share the same vocabulary in the same order, and their
#' document sets must be disjoint; documents are identified by row names
#' (or `quanteda::docid()` for dfms). In-memory shards are validated at
#' construction; path shards are validated when first materialized. The
#' functions that consume a corpus stream it one shard at a time, so peak
#' memory is one shard plus the per-document result vectors, whose length
#' is the total number of documents.
#'
#' Whether the corpus must hold counts or proportions is decided by the
#' consumer, exactly as for a single matrix: [optop_select()] expects
#' proportions, the partition, baseline, and index family expect counts.
#'
#' @return An object of class `optop_corpus`.
#'
#' @examples
#' # shard a small corpus in memory: results match the unsharded matrix
#' rdirich <- function(n, k) {
#'   g <- matrix(stats::rgamma(n * k, shape = 1), n, k)
#'   g / rowSums(g)
#' }
#' set.seed(1)
#' theta <- rdirich(40, 3)
#' rownames(theta) <- sprintf("d%02d", 1:40)
#' phi <- rdirich(3, 60)
#' colnames(phi) <- sprintf("w%02d", 1:60)
#' counts <- sim_dfm(theta, phi, doc_length = 400, seed = 2)
#' models <- list(optop_model(theta, phi))
#'
#' corp <- optop_corpus(list(counts[1:25, ], counts[26:40, ]))
#' corp
#'
#' part_sharded <- optop_make_partition(models, corp, c = 1)
#' part_plain <- optop_make_partition(models, counts, c = 1)
#' identical(part_sharded$nonrare_words, part_plain$nonrare_words)
#'
#' \dontrun{
#' # at real scale, point at per-shard files written upstream: only one
#' # shard is ever in memory
#' corp <- optop_corpus(list.files("shards", full.names = TRUE),
#'                      reader = function(p) qs2::qs_read(p))
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' @seealso [optop_select()], [optop_make_partition()],
#'   [optop_make_baseline()], [optop_index_holdout()]
#'
#' @export
optop_corpus <- function(x, reader = NULL) {
  if (inherits(x, "optop_corpus")) {
    return(x)
  }
  lazy <- is.character(x)
  if (lazy) {
    if (!length(x)) stop("x must name at least one shard file")
    if (!is.function(reader)) {
      stop("path shards need a reader function, e.g. readRDS")
    }
    missing_files <- !file.exists(x)
    if (any(missing_files)) {
      stop("shard file not found: ", x[which(missing_files)[1L]])
    }
    shards <- as.list(x)
  } else {
    shards <- if (is.list(x)) x else list(x)
    if (!length(shards)) stop("x must contain at least one shard")
    ok <- vapply(shards, function(s) {
      quanteda::is.dfm(s) || methods::is(s, "Matrix") || is.matrix(s)
    }, logical(1))
    if (!all(ok)) {
      stop("every shard must be a dfm or a (sparse) matrix")
    }
  }

  out <- structure(
    list(shards = shards, reader = reader, lazy = lazy,
         n_shards = length(shards), vocab = NULL),
    class = "optop_corpus"
  )

  if (!lazy) {
    vocab <- colnames(shards[[1L]])
    if (is.null(vocab)) {
      stop("shards must have column names (the vocabulary)")
    }
    same <- vapply(shards, function(s) identical(colnames(s), vocab),
                   logical(1))
    if (!all(same)) {
      stop("all shards must share the same vocabulary in the same order")
    }
    ids <- unlist(lapply(shards, .optop_shard_ids), use.names = FALSE)
    if (anyNA(ids) || anyDuplicated(ids)) {
      stop(paste("shard documents must carry unique identifiers (row names",
                 "or quanteda docids), disjoint across shards"))
    }
    out$vocab <- vocab
  }
  out
}

#' @export
print.optop_corpus <- function(x, ...) {
  kind <- if (x$lazy) "path" else "in-memory"
  cat(sprintf("An optop_corpus of %d %s shard%s\n", x$n_shards, kind,
              if (x$n_shards == 1L) "" else "s"))
  if (!is.null(x$vocab)) {
    cat(sprintf("  vocabulary: %d features\n", length(x$vocab)))
  }
  invisible(x)
}

# Document identifiers of one materialized shard
.optop_shard_ids <- function(s) {
  if (quanteda::is.dfm(s)) {
    as.character(quanteda::docid(s))
  } else {
    rownames(s)
  }
}

.optop_is_corpus <- function(x) inherits(x, "optop_corpus")

# Normalize any accepted dtm input to a corpus (a plain matrix becomes the
# trivial one-shard corpus, with no validation overhead beyond coercion)
.optop_as_corpus <- function(x, reader = NULL) {
  if (.optop_is_corpus(x)) x else optop_corpus(x, reader)
}

# Internal trivial wrap: a single matrix becomes a one-shard corpus without
# the constructor's naming requirements (legacy callers may pass unnamed
# dtms, which the single-shard code paths tolerate)
.optop_as_corpus_internal <- function(x) {
  if (.optop_is_corpus(x)) {
    return(x)
  }
  structure(list(shards = list(x), reader = NULL, lazy = FALSE,
                 n_shards = 1L, vocab = colnames(x)),
            class = "optop_corpus")
}

# Materialize shard i as a dgCMatrix with document identifiers as row names.
# Path shards are read through the corpus reader and validated here (the
# first materialization fixes the corpus vocabulary).
.optop_corpus_shard <- function(corpus, i) {
  s <- corpus$shards[[i]]
  if (corpus$lazy) {
    s <- corpus$reader(s)
    if (!(quanteda::is.dfm(s) || methods::is(s, "Matrix") || is.matrix(s))) {
      stop("the reader must return a dfm or a (sparse) matrix")
    }
  }
  ids <- .optop_shard_ids(s)
  m <- methods::as(s, "CsparseMatrix")
  if (!is.null(ids)) {
    rownames(m) <- ids
  }
  if (!is.null(corpus$vocab) && !identical(colnames(m), corpus$vocab)) {
    stop(sprintf("shard %d does not share the corpus vocabulary", i))
  }
  m
}

# Vocabulary of the corpus, materializing the first shard when lazy
.optop_corpus_vocab <- function(corpus) {
  if (!is.null(corpus$vocab)) {
    return(corpus$vocab)
  }
  colnames(.optop_corpus_shard(corpus, 1L))
}

# Materialize a fitted model supplied either directly or as a loader
# function (lazy models: the grid never needs to sit in memory at once)
.optop_materialize_model <- function(m) {
  if (is.function(m)) m() else m
}
