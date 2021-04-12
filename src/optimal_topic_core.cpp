#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

//' @keywords internal
// [[Rcpp::export]]
arma::mat optimal_topic_core(const Rcpp::List& lda_models,
                             const arma::sp_mat& weighted_dfm,
                             double q,
                             const Rcpp::CharacterVector& docs,
                             int n_docs,
                             int n_features)
{
    // this is the output of this method
    arma::mat Chi_K(lda_models.size(), 3);

    // loop LDA models
    arma::mat regstats(n_docs, 4);
    for (std::size_t i_mod = 0; i_mod < lda_models.size(); ++i_mod)
    {
        Rcpp::checkUserInterrupt();

        // getting the document word weights
        const Rcpp::S4& current_lda = lda_models[i_mod];
        const arma::mat& dtw = current_lda.slot("gamma");
        const arma::mat& beta = current_lda.slot("beta");
        arma::mat tww = exp(beta).t();
        int current_k = dtw.n_cols;
        // TODO TODO TODO
        // TODO current_k is useless, actually; there used to be a check at the
        // end when computing chi_out returns, that returned true for all rows;
        // I'm keeping current_k in Chi_J just because the calling function
        // expects a column "topic" in the matrix returned, and I don't know if
        // that's really needed and why
        // TODO TODO TODO

        Rcpp::Rcout << "---" << std::endl;
        Rcpp::Rcout << "# # # Processing LDA with k = " << current_k << std::endl;

        // looping over each document (k) in each model (j);
        // this is the loop that needs to be parallelized
        // Rcpp::Rcout << "--> Processing documents" << std::endl;
        for (std::size_t j_doc = 0; j_doc < n_docs; ++j_doc)
        {
            Rcpp::checkUserInterrupt();

            // subsetting dtw based on id_doc
            // TODO TODO TODO
            // TODO is there a way to make this sparse?? there should be; why are
            // we using all the zeros in the computation of chi_sq_fit below??
            // That way, we could use a reference to the row instead
            // TODO TODO TODO
            arma::vec weighted_dfm_j_doc(weighted_dfm.n_cols);
            const auto& row_j = weighted_dfm.row(j_doc);
            for (std::size_t row_index = 0; row_index < weighted_dfm_j_doc.size(); ++row_index)
            {
                weighted_dfm_j_doc(row_index) = row_j(row_index);
            }

            // element-wise multiplication
            // TODO this is the bottleneck now! But both tww and dtw are full, aren't they? There's
            // not much we can do
            arma::mat tww_dtw = tww.each_row() % dtw.row(j_doc);

            // get cumulative sum of sorted sums of rows of tww_dtw
            arma::vec X = arma::sum(tww_dtw, 1);
            // TODO second bottleneck (roughly half the time as tww_dtw)
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

            // compute the sum of the rest of the elements
            double sum_of_cut_weighted_dfm = arma::sum(weighted_dfm_j_doc.tail(n_cut_elements));
            double sum_of_cut_X = arma::sum(X.tail(n_cut_elements));
            double chi_sq_fit = icut * (
                arma::sum(
                          (
                           (weighted_dfm_j_doc.head(icut) - X.head(icut))
                           %
                           (weighted_dfm_j_doc.head(icut) - X.head(icut))
                          ) / X.head(icut)
                         )
                +
                (sum_of_cut_weighted_dfm - sum_of_cut_X) * (sum_of_cut_weighted_dfm - sum_of_cut_X) / sum_of_cut_X
            );

            regstats.row(j_doc) = arma::rowvec({double(current_k), double(j_doc), chi_sq_fit, double(icut)});
        }

        Chi_K.row(i_mod) = arma::rowvec({
            double(current_k),
            arma::sum(regstats.col(2)) / arma::sum(regstats.col(3)),
            R::pchisq(arma::sum(regstats.col(2)) / arma::sum(regstats.col(3)), 1, true, false)
        });
    }

    return Chi_K;
}
