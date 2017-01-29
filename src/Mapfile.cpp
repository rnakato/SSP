/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "Mapfile.hpp"

void SeqStatsGenome::readGenomeTable(const std::string &gt)
{
  std::vector<std::string> v;
  std::string lineStr;
  std::ifstream in(gt);
  if(!in) PRINTERR("Could nome open " << gt << ".");
  
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    SeqStats s(v[0], stoi(v[1]));
    chr.push_back(s);
  }
  return;
}

int32_t setIdLongestChr(SeqStatsGenome &genome)
{
  int32_t id(0);
  uint64_t lenmax(0);
  for(size_t i=0; i<genome.chr.size(); ++i) {
    if(lenmax < genome.chr[i].getlenmpbl()) {
      lenmax = genome.chr[i].getlenmpbl();
      id = i;
    }
  }
  return id;
}
