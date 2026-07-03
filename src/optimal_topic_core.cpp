#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

//' @keywords internal
// [[Rcpp::export]]
Rcpp::List optimal_topic_core(const Rcpp::List& lda_models,
                              const arma::sp_mat& weighted_dfm,
                              double q,
                              const arma::uvec& doc_map,
                              bool legacy,
                              bool return_envelope)
{
    const arma::uword n_models = lda_models.size();
    const arma::uword n_docs = weighted_dfm.n_rows;

    if (return_envelope && legacy) {
        Rcpp::stop("the envelope export is not available for the legacy pipeline");
    }
    if (return_envelope && n_models != 1) {
        Rcpp::stop("the envelope export works on a single model at a time");
    }

    // documents are processed in fixed-size blocks: one BLAS product per
    // (model, block) replaces the former dense W x K temporary per document,
    // while keeping memory bounded for large corpora
    const arma::uword block_size = 256;

    // transposed copy (W x J): a document becomes a contiguous CSC column, so
    // densifying a block of documents is a cheap column-range extraction
    // instead of element-by-element sparse row reads
    const arma::sp_mat dfm_t = weighted_dfm.t();

    // output: one row per model with the number of topics (the "topic" column
    // consumed by the caller), the raw Test 1 statistic of Eq. (8) in Lewis &
    // Grossetti (2022), and its degrees of freedom sum_j P_j; standardization
    // and p-values are computed by the R caller, where the tail convention is
    // visible and testable
    arma::mat Chi_K(n_models, 3);

    // per-document envelope, exported for the calibration layer: the bin
    // structure (kept-word probabilities plus the collapsed min-bin mass)
    // depends only on the model, so the null distribution of the statistic
    // can be simulated or moment-matched from these bins alone
    std::vector<double> bin_probs;
    std::vector<int> bin_counts;
    if (return_envelope) {
        bin_probs.reserve(n_docs * 128);
        bin_counts.reserve(n_docs);
    }

    const arma::uword max_gamma_row = doc_map.max();

    // loop LDA models
    for (arma::uword i_mod = 0; i_mod < n_models; ++i_mod)
    {
        Rcpp::checkUserInterrupt();

        // getting the document word weights
        const Rcpp::S4& current_lda = lda_models[i_mod];
        const arma::mat dtw = current_lda.slot("gamma");
        const arma::mat beta = current_lda.slot("beta");
        const arma::mat tww = arma::exp(beta).t();
        const arma::uword current_k = dtw.n_cols;

        if (max_gamma_row >= dtw.n_rows) {
            Rcpp::stop("document mapping points past the rows of @gamma: "
                       "the dfm and the LDA models are misaligned");
        }

        double sum_chi = 0.0;
        double sum_df = 0.0;

        for (arma::uword block_start = 0; block_start < n_docs; block_start += block_size)
        {
            Rcpp::checkUserInterrupt();

            const arma::uword block_end = std::min(block_start + block_size, n_docs) - 1;

            // W x b slabs, one column per document
            const arma::mat dfm_block(dfm_t.cols(block_start, block_end));
            const arma::uvec gamma_rows = doc_map.subvec(block_start, block_end);
            const arma::mat X_block = tww * arma::mat(dtw.rows(gamma_rows)).t();

            for (arma::uword j_col = 0; j_col < X_block.n_cols; ++j_col)
            {
                // get cumulative sum of sorted word probabilities
                arma::vec X = X_block.col(j_col);
                arma::vec weighted_dfm_j_doc = dfm_block.col(j_col);

                arma::uvec sort_indices = arma::sort_index(X, "descend");
                X = X(sort_indices);
                arma::vec X_cumsum = arma::cumsum(X);
                weighted_dfm_j_doc = weighted_dfm_j_doc(sort_indices);

                if (legacy) {
                    // 0.9.8 arithmetic, kept verbatim for reproducibility of
                    // pre-calibration results: rounded cutoff, crossing word
                    // collapsed into the tail, per-document multiplier P_j
                    auto first_greater_than_q = std::find_if(X_cumsum.begin(),
                                                             X_cumsum.end(),
                                                             [&q](double val) {return round(val * 1e4) / 1e4 > q;});
                    std::size_t icut = std::distance(X_cumsum.begin(), first_greater_than_q);
                    std::size_t n_cut_elements = X.size() - icut;

                    const arma::vec head_diff = weighted_dfm_j_doc.head(icut) - X.head(icut);
                    const double tail_X = arma::sum(X.tail(n_cut_elements));
                    const double tail_diff = arma::sum(weighted_dfm_j_doc.tail(n_cut_elements)) - tail_X;
                    const double pearson = arma::sum(head_diff % head_diff / X.head(icut))
                                           + tail_diff * tail_diff / tail_X;

                    sum_chi += double(icut) * pearson;
                    sum_df += double(icut);
                } else {
                    // Eq. (8): the P_j important words are the smallest head
                    // whose cumulative mass strictly exceeds q = 1 - I^K (the
                    // crossing word is kept, so the collapsed tail keeps mass
                    // strictly below I^K, cf. footnote 5), and the document
                    // contributes (bins) * Pearson with df = bins - 1, where
                    // bins = P_j + 1 with the collapsed min bin
                    auto first_greater_than_q = std::find_if(X_cumsum.begin(),
                                                             X_cumsum.end(),
                                                             [&q](double val) {return val > q;});
                    const std::size_t p_j = (first_greater_than_q == X_cumsum.end())
                        ? X.n_elem
                        : std::distance(X_cumsum.begin(), first_greater_than_q) + 1;
                    const std::size_t n_tail = X.n_elem - p_j;

                    const arma::vec head_diff = weighted_dfm_j_doc.head(p_j) - X.head(p_j);
                    double pearson = arma::sum(head_diff % head_diff / X.head(p_j));

                    std::size_t n_bins = p_j;
                    double tail_X = 0.0;
                    if (n_tail > 0) {
                        tail_X = arma::sum(X.tail(n_tail));
                        const double tail_diff = arma::sum(weighted_dfm_j_doc.tail(n_tail)) - tail_X;
                        pearson += tail_diff * tail_diff / tail_X;
                        n_bins += 1;
                    }

                    sum_chi += double(n_bins) * pearson;
                    sum_df += double(n_bins - 1);

                    if (return_envelope) {
                        for (std::size_t p = 0; p < p_j; ++p) {
                            bin_probs.push_back(X(p));
                        }
                        if (n_tail > 0) {
                            bin_probs.push_back(tail_X);
                        }
                        bin_counts.push_back(int(n_bins));
                    }
                }
            }
        }

        Chi_K.row(i_mod) = arma::rowvec({
            double(current_k),
            sum_chi,
            sum_df
        });
    }

    if (return_envelope) {
        return Rcpp::List::create(
            Rcpp::Named("stat") = Chi_K,
            Rcpp::Named("bin_probs") = Rcpp::NumericVector(bin_probs.begin(), bin_probs.end()),
            Rcpp::Named("bin_counts") = Rcpp::IntegerVector(bin_counts.begin(), bin_counts.end())
        );
    }
    return Rcpp::List::create(Rcpp::Named("stat") = Chi_K);
}
