/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "Mapfile.hpp"

void SeqStatsGenome::readGenomeTable(const std::string &gt, const int binsize) {
  std::vector<std::string> v;
  std::string lineStr;
  std::ifstream in(gt);
  if(!in) PRINTERR("Could nome open " << gt << ".");
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    SeqWigStats s(v[0], stoi(v[1]), binsize);
    chr.push_back(s);
  }
  return;
}
