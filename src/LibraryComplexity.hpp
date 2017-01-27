/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _LIBRARYCOMPLEXITY_HPP_
#define _LIBRARYCOMPLEXITY_HPP_

#include <iostream>
#include <string>
#include <stdint.h>
#include <unordered_map>
#include <boost/format.hpp>
#include "util.h"
#include "seq.h"
#include "SeqStats.hpp"
#include "BoostOptions.hpp"

class SeqStatsGenome;

class LibComp {
  bool pairedend;
  int32_t thredef;
  uint64_t ncmp;
  
  uint32_t nt_all, nt_nonred, nt_red;
  int32_t threshold;
  double r4cmp;
  bool lackOfRead;

 public:
 LibComp(const boost::program_options::variables_map &values):
  pairedend(values.count("pair")),
    thredef(values["thre_pb"].as<int32_t>()),
    ncmp(values["ncmp"].as<uint64_t>()),
    nt_all(0), nt_nonred(0), nt_red(0), threshold(0),
    r4cmp(0), lackOfRead(false) {}

  double getcomplexity() const { return getratio(nt_nonred, nt_all); }
  
  void setThreshold(const uint64_t nread, const uint64_t lenmpbl) {
    if(thredef) threshold = thredef;
    else {
      uint32_t thre = getratio(nread, lenmpbl) * 10;
      threshold = thre > 1 ? thre : 1;
    }
    std::cout << "Checking redundant reads: redundancy threshold " << threshold << std::endl;
  }
  int32_t getThreshold() const { return threshold; };

  void checkRedundantReads(SeqStatsGenome &genome);
  void filtering_eachchr_single(SeqStats &chr);
  void filtering_eachchr_pair(SeqStats &chr);
  void hashFilterCmpSingle(std::unordered_map<int32_t, int32_t> &mp, const SeqStats &chr, const Strand::Strand strand);
  void hashFilterCmpPair(std::unordered_map<std::string, int32_t> &mp, const SeqStats &chr, const Strand::Strand strand);

  void print(std::ofstream &out) const {
    if(lackOfRead) out << boost::format("Library complexity: (%1$.3f) (%2%/%3%)\n") % getcomplexity() % nt_nonred % nt_all;
    else           out << boost::format("Library complexity: %1$.3f (%2%/%3%)\n")   % getcomplexity() % nt_nonred % nt_all;
  }
};

#endif /* _LIBRARYCOMPLEXITY_HPP_ */
