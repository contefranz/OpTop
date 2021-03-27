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
    // Rcpp::Rcout << "ENTERED RCPP CODE" << std::endl;
  
    // this is the output of this method 
    arma::mat Chi_K(lda_models.size(), 3);

    // loop variables
    arma::mat regstats(n_docs, 4);
    for (std::size_t i_mod = 0; i_mod < lda_models.size(); ++i_mod)
    {
        Rcpp::checkUserInterrupt();

        // getting the document word weights
        Rcpp::S4 current_lda = lda_models[i_mod];
        arma::mat dtw = current_lda.slot("gamma");
        int current_k = dtw.n_cols;

        Rcpp::Rcout << "---" << std::endl;
        Rcpp::Rcout << "# # # Processing LDA with k = " << current_k << std::endl;

        // getting the term word weights -> beta
        arma::mat tww = exp(Rcpp::as<arma::mat>(current_lda.slot("beta"))).t();

        // adding row position to both objects
        // TODO could this be optimized further by creating these two matrices in
        // one go? or does armadillo take care of this seamlessly?
        arma::vec tmp1(dtw.n_rows);
        arma::vec tmp2(tww.n_rows);
        std::iota(tmp1.begin(), tmp1.end(), 1);
        std::iota(tmp2.begin(), tmp2.end(), 1);
        dtw.insert_cols(dtw.n_cols, tmp1);
        tww.insert_cols(tww.n_cols, tmp2);

        // looping over each document (k) in each model (j);
        // this is the loop that needs to be parallelized
        // Rcpp::Rcout << "--> Processing documents" << std::endl;
        for (std::size_t j_doc = 0; j_doc < n_docs; ++j_doc)
        {
            Rcpp::checkUserInterrupt();
            // TODO logic unclear; could all of this loop be written with fewer calls?

            // subsetting word proportions and dtw based on id_doc;
            // Rcpp::Rcout << "weighted_dfm row " << j_doc << ": " << std::endl << weighted_dfm.row(j_doc) << std::endl; // TODO REMOVE
            // define BestPair as n_cols(weighted_dfm) x 1, then fill it with weighted_dfm jth
            arma::mat BestPair(weighted_dfm.n_cols, 1);
            for (std::size_t row_index = 0; row_index < BestPair.n_rows; ++row_index)
            {
                BestPair(row_index, 0) = weighted_dfm.row(j_doc)(row_index);
            }
            // Rcpp::Rcout << "Best Pair: " << std::endl << BestPair << std::endl; // TODO REMOVE

            // casting N x K matrix
            arma::mat dtw_j_doc(n_features, dtw.n_cols);
            for (std::size_t row_index = 0; row_index < dtw_j_doc.n_rows; ++row_index)
            {
              dtw_j_doc.row(row_index) = dtw.row(j_doc);
            }
            // Rcpp::Rcout << "dtw_j_doc: " << std::endl << dtw_j_doc << std::endl; // TODO REMOVE

            // this avoids the use of j index which does not match with matlab code
            // in matlab j loops over k_start -> k_end
            // here starts from 1 up to the latest model
            arma::mat sub_dtw_j_doc = dtw_j_doc.cols(0, dtw_j_doc.n_cols - 2);
            // Rcpp::Rcout << "sub_dtw_j_doc: " << std::endl << sub_dtw_j_doc << std::endl; // TODO REMOVE
            arma::mat sub_tww = tww.cols(0, tww.n_cols - 2);
            // Rcpp::Rcout << "sub_tww: " << std::endl << sub_tww << std::endl; // TODO REMOVE

            // element-wise multiplication
            arma::mat tww_dtw = sub_dtw_j_doc % sub_tww;
            // Rcpp::Rcout << "tww_dtw: " << std::endl << tww_dtw << std::endl; // TODO REMOVE

            // this returns a vector...maybe we want a matrix
            arma::vec X = arma::sum(tww_dtw, 1);
            // Rcpp::Rcout << "X: " << std::endl << X << std::endl; // TODO REMOVE
            BestPair.insert_cols(BestPair.n_cols, X);
            BestPair = BestPair.rows(arma::sort_index(-BestPair.col(1)));
            // Rcpp::Rcout << "BestPair: " << std::endl << BestPair << std::endl; // TODO REMOVE

            // compute the cumulative probability over estimations
            BestPair.insert_cols(BestPair.n_cols, arma::cumsum(BestPair.col(1)));
            // Rcpp::Rcout << "BestPair: " << std::endl << BestPair << std::endl; // TODO REMOVE
            int n_BP = BestPair.n_rows;

            // stop when you reach q
            arma::mat AggBestPair = BestPair.rows(arma::find(arma::round(BestPair.col(2) * 1e4) / 1e4 <= q));
            // Rcpp::Rcout << "AggBestPair: " << std::endl << AggBestPair << std::endl; // TODO REMOVE
            int icut = AggBestPair.n_rows;
            AggBestPair.insert_rows(AggBestPair.n_rows, arma::sum(BestPair.rows(icut, n_BP - 1), 0));
            // Rcpp::Rcout << "AggBestPair: " << std::endl << AggBestPair << std::endl; // TODO REMOVE
            double chi_sq_fit = icut * arma::sum(
              (AggBestPair.col(0) - AggBestPair.col(1)) % (AggBestPair.col(0) - AggBestPair.col(1))
              / AggBestPair.col(1)
            );
            // Rcpp::Rcout << "chi_sq_fit: " << std::endl << chi_sq_fit << std::endl; // TODO REMOVE

            regstats.row(j_doc) = arma::rowvec({double(current_k), double(j_doc), chi_sq_fit, double(icut)});
            // Rcpp::Rcout << "regstats: " << std::endl << regstats << std::endl; // TODO REMOVE
        }

        arma::mat chi_out = regstats.rows(arma::find(regstats.col(0) == double(current_k)));
        // Rcpp::Rcout << "chi_out: " << std::endl << chi_out << std::endl; // TODO REMOVE
        Chi_K.row(i_mod) = arma::rowvec({
          double(current_k),
          arma::sum(chi_out.col(2)) / arma::sum(chi_out.col(3)),
          R::pchisq(arma::sum(chi_out.col(2)) / arma::sum(chi_out.col(3)), 1, true, false)
        });
        // Rcpp::Rcout << "Chi_K: " << std::endl << Chi_K << std::endl; // TODO REMOVE
    }

    // Rcpp::Rcout << "EXITING RCPP CODE" << std::endl;
    return Chi_K;
}