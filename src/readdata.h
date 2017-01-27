/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef READDATA_H
#define READDATA_H

#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "SeqStats.hpp"

using GeneDataMap = std::unordered_map<std::string, genedata>;
using HashOfGeneDataMap = std::unordered_map<std::string, GeneDataMap>;

int32_t countmp(HashOfGeneDataMap &);
std::vector<std::string> scanGeneName(const HashOfGeneDataMap &);
HashOfGeneDataMap extract_mp(const HashOfGeneDataMap &, const std::vector<std::string>);
std::vector<std::string> readGeneList(const std::string&);
HashOfGeneDataMap parseRefFlat(const std::string&);
HashOfGeneDataMap parseGtf(const std::string&);
HashOfGeneDataMap construct_gmp(const HashOfGeneDataMap &);
void printMap(const HashOfGeneDataMap &);
void printRefFlat(const HashOfGeneDataMap &, const int32_t nameflag);
std::vector<chrsize> read_genometable(const std::string&);

template <class T>
std::vector<T> parseBed(const std::string &fileName)
{
  std::vector<T> vbed;
  std::ifstream in(fileName);
  if(!in) PRINTERR("BED file does not exist.");

  std::string lineStr;
  std::vector<std::string> v;
  while (!in.eof()) {
    getline(in, lineStr);
    
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v[1] == "start") continue;
    T bed(v);
    vbed.push_back(bed);
  }

  return vbed;
}

template <class T>
void printBed(const std::vector<T> &vbed)
{
  for (auto &x: vbed) {
    x.print();
    std::cout << std::endl;
  }
  std::cout << "bed num: " << vbed.size() << std::endl;
  return;
}

#endif  // READDATA_H
