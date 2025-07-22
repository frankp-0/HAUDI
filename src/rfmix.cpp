#include <Rcpp.h>
#include <array>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

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

struct AncestryTract {
  std::string chrom;
  uint32_t spos;
  uint32_t epos;
  uint8_t ancestry0;
  uint8_t ancestry1;
};

// [[Rcpp::export]]
DataFrame rcpp_read_rfmix(std::string msp_file) {
  std::ifstream infile(msp_file);
  if (!infile.is_open()) {
    stop("Failed to open input file.");
  }

  std::string line;
  std::getline(infile, line); // skip population codes
  std::getline(infile, line); // header

  std::vector<std::string> hdr_fields;
  std::stringstream ss(line);
  std::string field;

  while (std::getline(ss, field, '\t')) {
    hdr_fields.push_back(field);
  }

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

  std::unordered_map<std::string, std::vector<AncestryTract>> sample_tracts;
  std::vector<uint8_t> prev_anc(n_hap, 0);
  std::vector<uint32_t> prev_spos(n_hap, 0);
  uint32_t cur_epos = 0;
  std::string cur_chrom = "chr0";
  bool is_first_crf = true;

  while (std::getline(infile, line)) {
    CRFPoint crf = CRFPoint::from_msp_line(line, n_hap);

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

    if (is_first_crf) {
      is_first_crf = false;
      prev_anc = crf.ancestries;
      prev_spos.assign(n_hap, crf.spos);
      cur_chrom = crf.chrom;
    } else {
      for (size_t i = 0; i < n_hap; i += 2) {
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

  for (size_t i = 0; i < n_hap; i += 2) {
    auto &[sample_name0, hap0] = hap_index_to_sample[i];
    auto &[sample_name1, hap1] = hap_index_to_sample[i + 1];
    assert(sample_name0 == sample_name1);

    sample_tracts[sample_name0].push_back(
        {cur_chrom, prev_spos[i], cur_epos, prev_anc[i], prev_anc[i + 1]});
  }

  // Flatten to DataFrame
  std::vector<std::string> samples, chroms;
  std::vector<uint32_t> spos_vec, epos_vec;
  std::vector<int> ancestry0_vec, ancestry1_vec;

  for (const auto &[sample, tracts] : sample_tracts) {
    for (const auto &tract : tracts) {
      samples.push_back(sample);
      chroms.push_back(tract.chrom);
      spos_vec.push_back(tract.spos);
      epos_vec.push_back(tract.epos);
      ancestry0_vec.push_back(tract.ancestry0);
      ancestry1_vec.push_back(tract.ancestry1);
    }
  }

  return DataFrame::create(_["sample"] = samples,
                           _["chrom"] = chroms,
                           _["spos"] = spos_vec,
                           _["epos"] = epos_vec,
                           _["ancestry0"] = ancestry0_vec,
                           _["ancestry1"] = ancestry1_vec);
}

