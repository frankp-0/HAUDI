#include <Rcpp.h>
#include <array>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// Struct to represent a CRF point in RFMix msp.tsv file
struct CRFPoint {
  std::string chrom;
  uint32_t spos;
  uint32_t epos;
  std::vector<uint8_t> ancestries;

  static CRFPoint from_msp_line(const std::string &line, size_t n_hap) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;

    while (std::getline(ss, field, '\t')) {
      fields.push_back(field);
    }

    assert(fields.size() >= 7);

    CRFPoint crf;
    crf.chrom = fields[0];
    crf.spos = std::stoul(fields[1]);
    crf.epos = std::stoul(fields[2]);

    for (size_t i = 6; i < fields.size(); ++i) {
      crf.ancestries.push_back(static_cast<uint8_t>(std::stoi(fields[i])));
    }

    assert(crf.ancestries.size() == n_hap);
    return crf;
  }
};

// Struct to represent an ancestry tract
struct AncestryTract {
  std::string chrom;
  uint32_t spos;
  uint32_t epos;
  uint8_t anc0;
  uint8_t anc1;
};

// Constructs data frame with ancestry tracts from RFMix input
//
// Parameters:
// - msp_file A string with the file path for an RFMix msp.tsv file
//
// Returns:
// A DataFrame where each row is an ancestry tract, with columns:
// - "sample": The sample ID
// - "chrom": The chromosome
// - "spos": The start position of the tract
// - "epos": The end position of the tract
// - "anc0": The ancestry at haplotype 0
// - "anc1": The ancestry at haplotype 1
// [[Rcpp::export]]
DataFrame rcpp_read_rfmix(std::string msp_file) {
  std::ifstream infile(msp_file);
  if (!infile.is_open()) {
    stop("Failed to open input file.");
  }

  std::string line;
  std::getline(infile, line); // skip population codes
  std::getline(infile, line); // header

  // Get header fields
  std::vector<std::string> hdr_fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) {
    hdr_fields.push_back(field);
  }

  // Initialize object that can convert from column index
  // to sample/haplotype0
  std::vector<std::pair<std::string, size_t>> hap_index_to_sample;
  for (size_t i = 6; i < hdr_fields.size(); ++i) {
    auto &hap_field = hdr_fields[i];
    size_t dot_pos = hap_field.rfind('.');
    assert(dot_pos != std::string::npos);

    std::string sample = hap_field.substr(0, dot_pos);
    size_t hap_idx = std::stoul(hap_field.substr(dot_pos + 1));
    hap_index_to_sample.emplace_back(sample, hap_idx);
  }
  size_t n_hap = hap_index_to_sample.size();
  assert(n_hap % 2 == 0); // must be even for hap0 and hap1 pairs

  // Initialize ancestry tracts
  std::unordered_map<std::string, std::vector<AncestryTract>> sample_tracts;
  std::vector<uint8_t> prev_anc(n_hap, 0);
  std::vector<uint32_t> prev_spos(n_hap, 0);
  uint32_t cur_epos = 0;
  std::string cur_chrom = "chr0";
  bool is_first_crf = true;

  // Loop through lines in RFMix file
  while (std::getline(infile, line)) {
    // Construct CRF point
    CRFPoint crf = CRFPoint::from_msp_line(line, n_hap);

    // Flush tracts if new chromosome (this shouldn't happen)
    if (crf.chrom != cur_chrom && !is_first_crf) {
      for (size_t i = 0; i < n_hap; i += 2) {
        auto &[sample_name0, hap0] = hap_index_to_sample[i];
        auto &[sample_name1, hap1] = hap_index_to_sample[i + 1];
        assert(sample_name0 == sample_name1);

        sample_tracts[sample_name0].push_back(
            {cur_chrom, prev_spos[i], crf.spos - 1, prev_anc[i], prev_anc[i + 1]});
      }
      prev_spos.assign(n_hap, crf.spos);
      prev_anc = crf.ancestries;
      cur_chrom = crf.chrom;
    }

    cur_epos = crf.epos;

    // Start tracts if first CRF point
    if (is_first_crf) {
      is_first_crf = false;
      prev_anc = crf.ancestries;
      prev_spos.assign(n_hap, crf.spos);
      cur_chrom = crf.chrom;
    } else {
      for (size_t i = 0; i < n_hap; i += 2) {
        // Close out tract and start new one if ancestry switches
        // for either haplotype
        if (crf.ancestries[i] != prev_anc[i] ||
            crf.ancestries[i + 1] != prev_anc[i + 1]) {
          auto &[sample_name0, hap0] = hap_index_to_sample[i];
          auto &[sample_name1, hap1] = hap_index_to_sample[i + 1];
          assert(sample_name0 == sample_name1);

          sample_tracts[sample_name0].push_back(
              {cur_chrom, prev_spos[i], crf.spos - 1, prev_anc[i], prev_anc[i + 1]});

          prev_spos[i] = crf.spos;
          prev_spos[i + 1] = crf.spos;
          prev_anc[i] = crf.ancestries[i];
          prev_anc[i + 1] = crf.ancestries[i + 1];
        }
      }
    }
  }

  // Close out all open tracts after file ends
  for (size_t i = 0; i < n_hap; i += 2) {
    auto &[sample_name0, hap0] = hap_index_to_sample[i];
    auto &[sample_name1, hap1] = hap_index_to_sample[i + 1];
    assert(sample_name0 == sample_name1);

    sample_tracts[sample_name0].push_back(
        {cur_chrom, prev_spos[i], cur_epos, prev_anc[i], prev_anc[i + 1]});
  }

  // Return tracts in DataFrame
  std::vector<std::string> samples, chroms;
  std::vector<uint32_t> spos_vec, epos_vec;
  std::vector<int> anc0_vec, anc1_vec;

  for (const auto &[sample, tracts] : sample_tracts) {
    for (const auto &tract : tracts) {
      samples.push_back(sample);
      chroms.push_back(tract.chrom);
      spos_vec.push_back(tract.spos);
      epos_vec.push_back(tract.epos);
      anc0_vec.push_back(tract.anc0);
      anc1_vec.push_back(tract.anc1);
    }
  }

  return DataFrame::create(_["sample"] = samples,
                           _["chrom"] = chroms,
                           _["spos"] = spos_vec,
                           _["epos"] = epos_vec,
                           _["anc0"] = anc0_vec,
                           _["anc1"] = anc1_vec);
}

