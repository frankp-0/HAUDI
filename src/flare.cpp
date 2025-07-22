#include <Rcpp.h>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>
using namespace Rcpp;

// ========================== Data Structures ==========================

struct AncestryTract {
  std::string chrom;
  uint32_t spos;
  uint32_t epos;
  uint8_t ancestry0;
  uint8_t ancestry1;
};

struct VCFRecord {
  std::string chrom;
  uint32_t pos;
  std::vector<uint8_t> ancestry_hap0;
  std::vector<uint8_t> ancestry_hap1;

  VCFRecord(size_t n_samples)
      : ancestry_hap0(n_samples, 255), ancestry_hap1(n_samples, 255) {}
};

// reads line from FLARE VCF into buffer, adding it to line until line break
std::string gz_readline(gzFile file) {
  const size_t chunk_size = 65536;
  char buffer[chunk_size];
  std::string line;

  while (gzgets(file, buffer, chunk_size)) {
    line += buffer;
    if (!line.empty() && line.back() == '\n')
      break;
  }

  while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
    line.pop_back();
  }
  return line;
}

// split line into fields/values
std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

// extract AN1 and AN2 for a sample
std::pair<uint8_t, uint8_t>
extract_AN1_AN2(const std::string &sample_field,
                const std::vector<std::string> &format_fields, int an1_idx,
                int an2_idx) {
  uint8_t missing = 255;
  std::vector<std::string> tokens = split(sample_field, ':');
  uint8_t an1 = missing, an2 = missing;

  if (an1_idx >= 0 && an1_idx < (int)tokens.size()) {
    try {
      int val = std::stoi(tokens[an1_idx]);
      if (val >= 0 && val <= 255)
        an1 = val;
    } catch (...) {}
  }
  if (an2_idx >= 0 && an2_idx < (int)tokens.size()) {
    try {
      int val = std::stoi(tokens[an2_idx]);
      if (val >= 0 && val <= 255)
        an2 = val;
    } catch (...) {}
  }
  return {an1, an2};
}

// build VCF record (chrom, pos, ancestries)
VCFRecord construct_vcf_record(const std::vector<std::string> &fields,
                               const std::vector<std::string> &format_fields,
                               int an1_idx, int an2_idx, int format_idx,
                               int n_samples, const std::string &chrom,
                               uint32_t pos) {
  VCFRecord record(n_samples);
  record.chrom = chrom;
  record.pos = pos;
  for (int i = 0; i < n_samples; ++i) {
    auto [an1, an2] = extract_AN1_AN2(fields[format_idx + 1 + i], format_fields,
                                      an1_idx, an2_idx);
    record.ancestry_hap0[i] = an1;
    record.ancestry_hap1[i] = an2;
  }
  return record;
}

// close all open tracts (called at end of chromosome)
void finalize_open_tracts(
    const std::vector<std::string> &sample_ids,
    const std::vector<uint8_t> &prev_anc,
    const std::vector<uint32_t> &prev_spos,
    std::unordered_map<std::string, std::vector<AncestryTract>> &sample_tracts,
    const std::string &chrom, uint32_t final_pos) {
  for (size_t i = 0; i < sample_ids.size(); ++i) {
    const std::string &sample = sample_ids[i];
    uint8_t anc0 = prev_anc[i * 2];
    uint8_t anc1 = prev_anc[i * 2 + 1];
    if (anc0 != 255 || anc1 != 255) {
      sample_tracts[sample].push_back(
        {chrom, prev_spos[i * 2], final_pos, anc0, anc1});
    }
  }
}

// Constructs data frame with ancestry tracts from FLARE input
// [[Rcpp::export]]
DataFrame rcpp_read_flare(std::string flare_file) {
  gzFile file = gzopen(flare_file.c_str(), "rb");
  if (!file) {
    std::cerr << "Failed to open input VCF file.\n";
    return DataFrame::create();  // Return empty DataFrame
  }

  std::string line;
  bool found_header = false;
  while (!(line = gz_readline(file)).empty()) {
    if (line.substr(0, 6) == "#CHROM") {
      found_header = true;
      break;
    }
  }
  if (!found_header) {
    std::cerr << "Missing #CHROM header line.\n";
    gzclose(file);
    return DataFrame::create();
  }

  std::vector<std::string> header_fields = split(line, '\t');
  int chrom_idx = -1, pos_idx = -1, format_idx = -1;
  for (size_t i = 0; i < header_fields.size(); ++i) {
    if (header_fields[i] == "#CHROM")
      chrom_idx = i;
    else if (header_fields[i] == "POS")
      pos_idx = i;
    else if (header_fields[i] == "FORMAT")
      format_idx = i;
  }
  if (chrom_idx == -1 || pos_idx == -1 || format_idx == -1) {
    std::cerr << "Missing essential VCF columns.\n";
    gzclose(file);
    return DataFrame::create();
  }

  std::vector<std::string> sample_ids;
  for (size_t i = format_idx + 1; i < header_fields.size(); ++i) {
    sample_ids.push_back(header_fields[i]);
  }

  int n_samples = sample_ids.size();
  std::unordered_map<std::string, std::vector<AncestryTract>> sample_tracts;
  std::vector<uint8_t> prev_anc(n_samples * 2, 255);
  std::vector<uint32_t> prev_spos(n_samples * 2, 0);
  uint32_t cur_pos = 0;
  uint32_t prev_pos = 0;
  std::string cur_chrom = "chr0";
  bool is_first_record = true;

  while (!(line = gz_readline(file)).empty()) {
    if (line.empty() || line[0] == '#')
      continue;

    std::vector<std::string> fields = split(line, '\t');
    if (fields.size() < format_idx + 1 + n_samples)
      continue;

    std::string chrom = fields[chrom_idx];
    uint32_t pos = std::stoul(fields[pos_idx]);

    std::vector<std::string> format_fields = split(fields[format_idx], ':');
    int an1_idx = -1, an2_idx = -1;
    for (size_t i = 0; i < format_fields.size(); ++i) {
      if (format_fields[i] == "AN1")
        an1_idx = i;
      else if (format_fields[i] == "AN2")
        an2_idx = i;
    }

    VCFRecord record =
        construct_vcf_record(fields, format_fields, an1_idx, an2_idx,
                             format_idx, n_samples, chrom, pos);

    prev_pos = cur_pos;
    cur_pos = pos;

    if (is_first_record) {
      for (size_t i = 0; i < n_samples; ++i) {
        prev_spos[i * 2] = pos;
        prev_spos[i * 2 + 1] = pos;
        prev_anc[i * 2] = record.ancestry_hap0[i];
        prev_anc[i * 2 + 1] = record.ancestry_hap1[i];
      }
      cur_chrom = chrom;
      is_first_record = false;
      continue;
    }

    if (chrom != cur_chrom) {
      finalize_open_tracts(sample_ids, prev_anc, prev_spos, sample_tracts,
                           cur_chrom, pos - 1);
      for (size_t i = 0; i < n_samples; ++i) {
        prev_spos[i * 2] = pos;
        prev_spos[i * 2 + 1] = pos;
        prev_anc[i * 2] = record.ancestry_hap0[i];
        prev_anc[i * 2 + 1] = record.ancestry_hap1[i];
      }
      cur_chrom = chrom;
      continue;
    }

    for (size_t i = 0; i < n_samples; ++i) {
      const std::string &sample = sample_ids[i];
      size_t idx0 = i * 2;
      size_t idx1 = i * 2 + 1;
      uint8_t new_anc0 = record.ancestry_hap0[i];
      uint8_t new_anc1 = record.ancestry_hap1[i];
      if ((new_anc0 != prev_anc[idx0]) || (new_anc1 != prev_anc[idx1])) {
        uint32_t midpoint = prev_pos + (cur_pos - prev_pos) / 2;
        sample_tracts[sample].push_back({
          chrom, prev_spos[idx0], midpoint, prev_anc[idx0], prev_anc[idx1]
        });
        prev_spos[idx0] = midpoint + 1;
        prev_anc[idx0] = new_anc0;
        prev_anc[idx1] = new_anc1;
      }
    }
  }
  gzclose(file);

  finalize_open_tracts(sample_ids, prev_anc, prev_spos, sample_tracts,
                       cur_chrom, cur_pos);

  std::vector<std::string> samples, chroms;
  std::vector<uint32_t> spos_vec, epos_vec;
  std::vector<int> ancestry_hap0_vec, ancestry_hap1_vec;

  for (const auto &[sample, tracts] : sample_tracts) {
    for (const auto &tract : tracts) {
      samples.push_back(sample);
      chroms.push_back(tract.chrom);
      spos_vec.push_back(tract.spos);
      epos_vec.push_back(tract.epos);
      ancestry_hap0_vec.push_back(static_cast<int>(tract.ancestry0));
      ancestry_hap1_vec.push_back(static_cast<int>(tract.ancestry1));
    }
  }

  return DataFrame::create(
      _["sample"] = samples,
      _["chrom"] = chroms,
      _["spos"] = spos_vec,
      _["epos"] = epos_vec,
      _["ancestry0"] = ancestry_hap0_vec,
      _["ancestry1"] = ancestry_hap1_vec
  );
}
