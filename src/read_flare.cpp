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
  uint8_t ancestry;
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

  // keep adding buffer to new line until end of FLARE line is reached
  while (gzgets(file, buffer, chunk_size)) {
    line += buffer;
    if (!line.empty() && line.back() == '\n')
      break;
  }

  // strip \n from line
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
    } catch (...) {
    }
  }
  if (an2_idx >= 0 && an2_idx < (int)tokens.size()) {
    try {
      int val = std::stoi(tokens[an2_idx]);
      if (val >= 0 && val <= 255)
        an2 = val;
    } catch (...) {
    }
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
    std::unordered_map<std::string, std::array<std::vector<AncestryTract>, 2>>
        &sample_tracts,
    const std::string &chrom, uint32_t final_pos) {
  for (size_t i = 0; i < sample_ids.size(); ++i) {
    const std::string &sample = sample_ids[i];
    for (size_t hap = 0; hap < 2; ++hap) {
      size_t idx = i * 2 + hap;
      if (prev_anc[idx] != 255) {
        sample_tracts[sample][hap].push_back(
            {chrom, prev_spos[idx], final_pos, prev_anc[idx]});
      }
    }
  }
}

// Constructs data frame with ancestry tracts from FLARE input
// [[Rcpp::export]]
DataFrame rcpp_read_flare(std::string flare_file) {
  gzFile file = gzopen(flare_file.c_str(), "rb");
  if (!file) {
    std::cerr << "Failed to open input VCF file.\n";
    return 1;
  }

  // skip header lines until column definition
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
    return 1;
  }

  // process header, extracting index of chromosome,
  // position, and format on future lines
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
    return 1;
  }

  // extract sample IDs
  std::vector<std::string> sample_ids;
  for (size_t i = format_idx + 1; i < header_fields.size(); ++i) {
    sample_ids.push_back(header_fields[i]);
  }

  // initialize values
  int n_samples = sample_ids.size();
  std::unordered_map<std::string, std::array<std::vector<AncestryTract>, 2>>
      sample_tracts;
  std::vector<uint8_t> prev_anc(n_samples * 2, 255);
  std::vector<uint32_t> prev_spos(n_samples * 2, 0);
  uint32_t cur_pos = 0;
  uint32_t prev_pos = 0;
  std::string cur_chrom = "chr0";
  bool is_first_record = true;

  while (!(line = gz_readline(file)).empty()) {
    if (line.empty() || line[0] == '#')
      continue;

    // split line into elements
    std::vector<std::string> fields = split(line, '\t');
    if (fields.size() < format_idx + 1 + n_samples)
      continue;

    // get current chromosome, position
    std::string chrom = fields[chrom_idx];
    uint32_t pos = std::stoul(fields[pos_idx]);

    std::vector<std::string> format_fields = split(fields[format_idx], ':');
    // check where ancestry fields occur
    int an1_idx = -1, an2_idx = -1;
    for (size_t i = 0; i < format_fields.size(); ++i) {
      if (format_fields[i] == "AN1")
        an1_idx = i;
      else if (format_fields[i] == "AN2")
        an2_idx = i;
    }

    // build current VCF record
    VCFRecord record =
        construct_vcf_record(fields, format_fields, an1_idx, an2_idx,
                             format_idx, n_samples, chrom, pos);

    // record previous and current record's position
    prev_pos = cur_pos;
    cur_pos = pos;

    // update info and start new tracts if first record
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

    // close open tracts and update info if new chromosome
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

    // push new tracts and update spos, prev_anc where ancestry changes
    for (size_t i = 0; i < n_samples; ++i) {
      const std::string &sample = sample_ids[i];
      for (size_t hap = 0; hap < 2; ++hap) {
        size_t idx = i * 2 + hap;
        uint8_t new_anc =
            hap == 0 ? record.ancestry_hap0[i] : record.ancestry_hap1[i];
        if (new_anc != prev_anc[idx]) {
          uint32_t midpoint = prev_pos + (cur_pos - prev_pos) / 2;
          sample_tracts[sample][hap].push_back(
              {chrom, prev_spos[idx], midpoint, prev_anc[idx]});
          prev_spos[idx] = midpoint + 1;
          prev_anc[idx] = new_anc;
        }
      }
    }
  }
  gzclose(file);

  // close out tracts when file ends
  finalize_open_tracts(sample_ids, prev_anc, prev_spos, sample_tracts,
                       cur_chrom, cur_pos);

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
