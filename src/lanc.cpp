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

// [[Rcpp::export]]
NumericMatrix rcpp_get_masked_geno(IntegerMatrix anc0, IntegerMatrix anc1,
                                   NumericMatrix gen0, NumericMatrix gen1,
                                   int n_anc) {
  int n_samp = anc0.nrow();
  int n_var = anc0.ncol();

  NumericMatrix out(n_samp, n_var * n_anc + n_var);

  const int *p_anc0 = INTEGER(anc0);
  const int *p_anc1 = INTEGER(anc1);
  const double *p_gen0 = REAL(gen0);
  const double *p_gen1 = REAL(gen1);
  double *p_out = REAL(out);

  // ancestry blocks
  for (int anc = 0; anc < n_anc; ++anc) {
    int offset = anc * n_var;

    for (int j = 0; j < n_var; ++j) {

      const int *h0 = p_anc0 + j * n_samp;
      const int *h1 = p_anc1 + j * n_samp;
      const double *g0 = p_gen0 + j * n_samp;
      const double *g1 = p_gen1 + j * n_samp;

      double *out_col = p_out + (offset + j) * n_samp;

      for (int i = 0; i < n_samp; ++i) {
        double v = 0.0;
        if (h0[i] == anc)
          v += g0[i];
        if (h1[i] == anc)
          v += g1[i];
        out_col[i] = v;
      }
    }
  }

  int final_offset = n_var * n_anc;

  for (int j = 0; j < n_var; ++j) {
    const double *g0 = p_gen0 + j * n_samp;
    const double *g1 = p_gen1 + j * n_samp;
    double *out_col = p_out + (final_offset + j) * n_samp;

    for (int i = 0; i < n_samp; ++i)
      out_col[i] = g0[i] + g1[i];
  }

  return out;
}
