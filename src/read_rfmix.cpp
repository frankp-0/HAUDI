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
  uint8_t ancestry;
};

// [[Rcpp::export]]
DataFrame rcpp_read_rfmix(std::string msp_file) {
  std::ifstream infile(msp_file);
  if (!infile.is_open()) {
    std::cerr << "Failed to open input file.\n";
    return 1;
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

  std::unordered_map<std::string, std::array<std::vector<AncestryTract>, 2>>
      sample_tracts;
  std::vector<uint8_t> prev_anc(n_hap, 0);
  std::vector<uint32_t> prev_spos(n_hap, 0);
  uint32_t cur_epos = 0;
  std::string cur_chrom = "chr0";
  bool is_first_crf = true;

  while (std::getline(infile, line)) {
    CRFPoint crf = CRFPoint::from_msp_line(line, n_hap);

    if (crf.chrom != cur_chrom && !is_first_crf) {
      for (size_t i = 0; i < n_hap; ++i) {
        auto &[sample_name, hap_idx] = hap_index_to_sample[i];
        sample_tracts[sample_name][hap_idx].push_back(
            {cur_chrom, prev_spos[i], crf.spos - 1, prev_anc[i]});
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
      for (size_t i = 0; i < n_hap; ++i) {
        if (crf.ancestries[i] != prev_anc[i]) {
          auto &[sample_name, hap_idx] = hap_index_to_sample[i];
          sample_tracts[sample_name][hap_idx].push_back(
              {cur_chrom, prev_spos[i], crf.spos - 1, prev_anc[i]});
          prev_spos[i] = crf.spos;
          prev_anc[i] = crf.ancestries[i];
        }
      }
    }
  }

  for (size_t i = 0; i < n_hap; ++i) {
    auto &[sample_name, hap_idx] = hap_index_to_sample[i];
    sample_tracts[sample_name][hap_idx].push_back(
        {cur_chrom, prev_spos[i], cur_epos, prev_anc[i]});
  }

  // Flatten and return as DataFrame
  std::vector<std::string> samples, chroms;
  std::vector<int> haps;
  std::vector<uint32_t> spos_vec, epos_vec;
  std::vector<int> ancestry_vec;

  for (const auto &[sample, tracts_by_hap] : sample_tracts) {
    for (size_t hap = 0; hap < 2; ++hap) {
      for (const auto &tract : tracts_by_hap[hap]) {
        samples.push_back(sample);
        haps.push_back(hap);
        chroms.push_back(tract.chrom);
        spos_vec.push_back(tract.spos);
        epos_vec.push_back(tract.epos);
        ancestry_vec.push_back(tract.ancestry);
      }
    }
  }

  return DataFrame::create(_["sample"] = samples, _["hap"] = haps,
                           _["chrom"] = chroms, _["spos"] = spos_vec,
                           _["epos"] = epos_vec, _["ancestry"] = ancestry_vec);
}
