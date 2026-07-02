#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

//' @keywords internal
// [[Rcpp::export]]
arma::mat optimal_topic_core(const Rcpp::List& lda_models,
                             const arma::sp_mat& weighted_dfm,
                             double q,
                             const arma::uvec& doc_map)
{
    const arma::uword n_models = lda_models.size();
    const arma::uword n_docs = weighted_dfm.n_rows;

    // documents are processed in fixed-size blocks: one BLAS product per
    // (model, block) replaces the former dense W x K temporary per document,
    // while keeping memory bounded for large corpora
    const arma::uword block_size = 256;

    // transposed copy (W x J): a document becomes a contiguous CSC column, so
    // densifying a block of documents is a cheap column-range extraction
    // instead of element-by-element sparse row reads
    const arma::sp_mat dfm_t = weighted_dfm.t();

    // this is the output of this method: one row per model with the number of
    // topics (the "topic" column consumed by the caller), the standardized
    // chi-square statistic and its p-value
    arma::mat Chi_K(n_models, 3);

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
        double sum_icut = 0.0;

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

                // find elements with estimated probability lower than quantile
                auto first_greater_than_q = std::find_if(X_cumsum.begin(),
                                                         X_cumsum.end(),
                                                         [&q](double val) {return round(val * 1e4) / 1e4 > q;});
                std::size_t icut = std::distance(X_cumsum.begin(), first_greater_than_q);
                std::size_t n_cut_elements = X.size() - icut;

                // Pearson terms on the kept elements plus one pooled term for
                // the truncated tail
                const arma::vec head_diff = weighted_dfm_j_doc.head(icut) - X.head(icut);
                const double tail_X = arma::sum(X.tail(n_cut_elements));
                const double tail_diff = arma::sum(weighted_dfm_j_doc.tail(n_cut_elements)) - tail_X;
                const double chi_sq_fit = icut * (
                    arma::sum(head_diff % head_diff / X.head(icut))
                    +
                    tail_diff * tail_diff / tail_X
                );

                sum_chi += chi_sq_fit;
                sum_icut += double(icut);
            }
        }

        const double stat = sum_chi / sum_icut;
        Chi_K.row(i_mod) = arma::rowvec({
            double(current_k),
            stat,
            R::pchisq(stat, 1, true, false)
        });
    }

    return Chi_K;
}
