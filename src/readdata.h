/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef READGENE_H
#define READGENE_H

#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "seq.h"
#include "macro.h"

using GeneDataMap = std::unordered_map<std::string, genedata>;
using HashOfGeneDataMap = std::unordered_map<std::string, GeneDataMap>;

enum bpstatus {UNMAPPABLE, INBED, MAPPABLE, COVREAD_ALL, COVREAD_NORM};

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

std::vector<int32_t> readMpbl(std::string, std::string, int32_t, int32_t);
std::vector<int8_t> readMpbl_binary(int32_t);
std::vector<int8_t> readMpbl_binary(std::string, std::string, int32_t);
std::vector<int8_t> arraySetBed(std::vector<int8_t> &, std::string, const std::vector<bed> &);
std::string IntToString(int64_t n);

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
  for (auto x: vbed) {
    x.print();
    std::cout << std::endl;
  }
  std::cout << "bed num: " << vbed.size() << std::endl;
  return;
}

#endif  // READGENE_H
