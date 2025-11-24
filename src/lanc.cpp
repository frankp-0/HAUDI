#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_query_tracts(const IntegerVector &left_haps,
                       const IntegerVector &right_haps,
                       const IntegerVector &breakpoints,
                       const IntegerVector &offsets,
                       const IntegerVector &indices) {
  int n_samples = offsets.size() - 1;
  int n_variants = indices.size();

  // Output: uint8 equivalents in R = raw vectors, but matrix of ints is fine
  // unless huge
  IntegerMatrix left_out(n_samples, n_variants);
  IntegerMatrix right_out(n_samples, n_variants);

  for (int i = 0; i < n_samples; i++) {
    int start = offsets[i];
    int end = offsets[i + 1];

    int j = 0;
    int end_len = end - start;

    for (int q = 0; q < n_variants; q++) {
      int idx = indices[q];

      // advance j while breakpoint <= idx
      while (j < end_len && idx >= breakpoints[start + j]) {
        j++;
      }

      // j selects the ancestry segment
      left_out(i, q) = left_haps[start + j];
      right_out(i, q) = right_haps[start + j];
    }
  }

  return List::create(_["hap0"] = left_out, _["hap1"] = right_out);
}
