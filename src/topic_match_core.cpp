#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <memory>
#include "utils.h"
// [Rcpp::depends(RcppArmadillo)]

// ---------------- AUXILIARY FUNCTIONS ---------------- //
// get the maximum size in a vector of vectors
std::size_t get_max_size(const std::vector<std::shared_ptr<std::vector<std::size_t>>>& v)
{
    std::size_t max_size = 0;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        max_size = std::max<std::size_t>(max_size, v[i]->size());
    }

    return max_size;
}


// build a matrix from data stored in a vector of vectors; the matrix is built
// row-wise; for each row, missing data are filled with NA;
// type changes from arma::uword to arma::sword (i.e. an arma::imat is returned
// instead of an arma::umat, even though the input vector is on unsigned
// integers), because the return matrix needs to store NA
arma::imat list_of_vectors_to_matrix(const std::vector<std::shared_ptr<std::vector<std::size_t>>>& v,
                                     std::size_t n_cols)
{
    std::size_t n_rows = v.size();
    arma::imat out(n_rows, n_cols);
    for (std::size_t i = 0; i < n_rows; ++i)
    {
        for (std::size_t j = 0; j < n_cols; ++j)
        {
            if (j < v[i]->size())
            {
                out(i, j) = (*(v[i]))[j];
            }
            else
            {
                out(i, j) = NA_INTEGER;
            }
        }
    }

    return out;
}


//' @keywords internal
// [[Rcpp::export]]
Rcpp::List topic_match_core(const Rcpp::List& lda_models,
                            std::size_t best_pos, // this is supposed to be in 1-indexing format, so the same as R
                            std::size_t optimal_k,
                            bool var_correction)
{
    // TODO R code used this check, but it makes no sense; optimal_model is
    // numeric, since it is checked above
    // if (!is.LDA_VEM (optimal_model))
    const Rcpp::S4& lda_best = lda_models[best_pos - 1]; // use -1 because best_pos is in 1-indexing format
    const arma::mat& beta_best = lda_best.slot("beta");
    arma::mat tww_best = exp(beta_best).t();
    arma::mat tww_best_normalized = normalize_columns(tww_best);

    std::size_t start_index = best_pos;
    std::size_t end_index = lda_models.size();

    std::vector<std::shared_ptr<std::vector<std::size_t>>> BestMatch;
    std::vector<std::shared_ptr<std::vector<std::size_t>>> LeastMatch;

    // reserve memory for BestMatch and LeastMatch
    BestMatch.reserve(end_index - start_index);
    LeastMatch.reserve(end_index - start_index);

    for (std::size_t i_mod = start_index; i_mod < end_index; ++i_mod)
    {
        Rcpp::checkUserInterrupt();

        // reading tww
        const Rcpp::S4& current_lda = lda_models[i_mod];
        const arma::mat& beta = current_lda.slot("beta");
        arma::mat tww = exp(beta).t();
        int current_k = tww.n_cols;

        Rcpp::Rcout << "---" << std::endl;
        Rcpp::Rcout << "# # # Processing LDA with k = " << current_k << std::endl;

        // normalizing by scaling by vector norms
        arma::mat tww_normalized = normalize_columns(tww);

        // this is the matrix multiplication to get the Cosine Similarities
        // between the best set of topics and the target ones
        arma::mat CosSim = tww_best_normalized.t() * tww_normalized;

        // make sure CosSim has no NaNs
        if (CosSim.has_nan())
        {
            Rcpp::stop("CosSim matrix is not supposed to contain NaNs. Aborting");
        }

        // get the indices of the maximum similarity factors, then sort it;
        // sorting will come in handy when the list of informative and uninformative
        // factors is created later on
        arma::uvec maxsim(CosSim.n_rows);
        for (std::size_t i = 0; i < CosSim.n_rows; ++i)
        {
            maxsim[i] = CosSim.row(i).index_max();
        }
        maxsim = arma::sort(maxsim);


        // Identify minimum cosine similarity threshold.
        // If CosSim exceeds thresh it is too similar to be treated as a
        // distinct factor. Matlab uses a weighting of 0.7 which is difficult
        // to replicate here;
        const arma::vec& CosSimVec = CosSim.as_col();
        double thresh;
        if (var_correction)
        {
            thresh = arma::mean(CosSimVec) + 2.58 * arma::stddev(CosSimVec);
        }
        else
        {
            std::size_t n = CosSimVec.size();
            thresh =  arma::mean(CosSimVec)
                + 2.58 * arma::stddev(CosSimVec) * sqrt((n - 1) / n);
        }

        // in this next loop, we identify factors that are highly similar to at
        // least one base model factor, then merge with the factors considered
        // to be most similar and identify those that are considered
        // uninformative;

        // create vectors to hold informative and uninformative factors, with a
        // capacity given by current_k + 1 (+1 because we always want them to
        // contain current_k as first element) TODO CHECK
        auto informative = std::make_shared<std::vector<std::size_t>>();
        auto uninformative = std::make_shared<std::vector<std::size_t>>();
        informative->reserve(current_k + 1);
        uninformative->reserve(current_k + 1);
        informative->push_back(current_k);
        uninformative->push_back(current_k);

        arma::uvec::const_iterator current_maxsim = maxsim.begin();
        for (std::size_t j_factor = 0; j_factor < current_k; ++j_factor)
        {
            // loop through maxsim until the current value is greater than or
            // equal to j_factor; since maxsim was sorted, this allows us to
            // know whether j_factor was in maxsim or not without calling an
            // expensive find() function
            while (*current_maxsim < j_factor)
            {
                current_maxsim++;
            }
            bool j_factor_is_in_maxsim = (*current_maxsim == j_factor);

            // check if this column of CosSim has any element greater than threshold
            bool found = false;
            for (std::size_t i_elem = 0; i_elem < CosSim.n_rows && !found; ++i_elem)
            {
                if (CosSim(i_elem, j_factor) > thresh)
                {
                    found = true;
                }
            }

            // if found or if j_factor is in maxsim, then add j_factor to the
            // list of informative factors; else add it to the list of
            // uninformative factors;
            // NB: j_factor is added to "informative" and "uninformative"
            // *incremented by 1*, since these two vectors will then directly be
            // used to build the return matrices, which will then be used by R
            // (where 1-based indexing is used)
            if (found or j_factor_is_in_maxsim)
            {
                informative->push_back(j_factor + 1);
            }
            else
            {
                uninformative->push_back(j_factor + 1);
            }
        }

        // store data in output variables
        BestMatch.push_back(informative);
        if (uninformative->size() > 1)
        {
            LeastMatch.push_back(uninformative);
        }
    }

    // converting lists to matrix accounting for different vectors lengths
    std::size_t best_match_max_length = get_max_size(BestMatch);
    std::size_t least_match_max_length = get_max_size(LeastMatch);

    arma::imat BestMatchOut = list_of_vectors_to_matrix(BestMatch,
                                                        best_match_max_length);
    arma::imat LeastMatchOut = list_of_vectors_to_matrix(LeastMatch,
                                                         least_match_max_length);

    Rcpp::List out;
    out["BestMatch"] = BestMatchOut;
    out["LeastMatch"] = LeastMatchOut;
    return out;
}
