#include <Rcpp.h>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace Rcpp;

// Parse lines (excluding header) from a .lanc file into a list with
// a DataFrame for each sample/line. 
//
// Parameters:
// - lines: A character vector where each element is a line (excluding
// the header) from a .lanc file
//
// Returns:
// A list of DataFrames. Each element of the list corresponds
// to a sample/line. DataFrames represent ancestry tracts
// and contain columns:
// - "index": The index in the plink2 file where this tract ends
// - "anc0": The ancestry for haplotype 0 on this tract
// - "anc1": The ancestry for haplotype 1 on this tract
// [[Rcpp::export]]
List rcpp_parse_lanc(CharacterVector lines) {
    int n_samples = lines.size();
    List result(n_samples);

    for (int i = 0; i < n_samples; ++i) {
        std::istringstream iss(Rcpp::as<std::string>(lines[i]));
        std::string token;

        std::vector<int> indices;
        std::vector<int> anc0s;
        std::vector<int> anc1s;

        while (iss >> token) {
            auto colon_pos = token.find(':');
            int idx = std::stoi(token.substr(0, colon_pos));
            int anc0 = token[colon_pos + 1] - '0';
            int anc1 = token[colon_pos + 2] - '0';

            indices.push_back(idx);
            anc0s.push_back(anc0);
            anc1s.push_back(anc1);
        }

        result[i] = DataFrame::create(
            Named("index") = indices,
            Named("anc0") = anc0s,
            Named("anc1") = anc1s
        );
    }

    return result;
}

// Query indices from a plink2 fileset for ancestry
//
// Parameters:
// - query_indices: A vector of indices in the plink2 fileset
// - tract_data: A list with ancestry tracts (as returned by rcpp_parse_lanc)
//
// Returns:
// A list with two matrices:
// - "hap0": ancestry values for each sample (rows) and index (columns) on haplotype 0 
// - "hap1": ancestry values for each sample (rows) and index (columns) on haplotype 1 
// [[Rcpp::export]]
List rcpp_query_tracts(IntegerVector query_indices, List tract_data) {
    int n_samples = tract_data.size();
    int n_query = query_indices.size();

    IntegerMatrix hap0_matrix(n_samples, n_query);
    IntegerMatrix hap1_matrix(n_samples, n_query);

    const int* p_queries = INTEGER(query_indices);

    for (int i = 0; i < n_samples; ++i) {
        DataFrame df = as<DataFrame>(tract_data[i]);
        IntegerVector indices = df["index"];
        IntegerVector anc0s = df["anc0"];
        IntegerVector anc1s = df["anc1"];

        const int* p_indices = INTEGER(indices);
        const int* p_anc0s = INTEGER(anc0s);
        const int* p_anc1s = INTEGER(anc1s);
        int n_tract = indices.size();

        for (int j = 0; j < n_query; ++j) {
            int q = p_queries[j];
            // Use std::lower_bound instead of forward linear pass
            // because there are few tracts per sample
            const int* it = std::lower_bound(p_indices, p_indices + n_tract, q);
            int idx = it - p_indices;
            hap0_matrix(i, j) = p_anc0s[idx];
            hap1_matrix(i, j) = p_anc1s[idx];
        }
    }

    return List::create(
        Named("hap0") = hap0_matrix,
        Named("hap1") = hap1_matrix
    );
}

