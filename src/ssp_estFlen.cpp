/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#include "ssp_estFlen.h"
#include "ssp_shiftprofile.h"
#include <time.h>

void estimateFragLength(const MyOpt::Variables &values, Mapfile &p)
{
  if(values.count("pair") || values.count("nomodel")) return;

  clock_t t1,t2;
  t1 = clock();
  strShiftProfile(values, p, "jaccard"); 
  t2 = clock();
  std::cout << "Jaccard Bit: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  makeFCSProfile(values, p, "fcs");
  p.printSSPstats();

  clock_t t3 = clock();
  std::cout << "Fragment variability: " << static_cast<double>(t3 - t2) / CLOCKS_PER_SEC << "sec.\n";

  if(values.count("ssp_exjac")) {
    t1 = clock();
    strShiftProfile(values, p, "exjaccard");
    t2 = clock();
    std::cout << "Jaccard Vec: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }

  if(values.count("ssp_hd")) {
    t1 = clock();
    strShiftProfile(values, p, "hdp");
    t2 = clock();
    std::cout << "Hamming: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }
    
  if(values.count("ssp_cc")) {
    t1 = clock();
    strShiftProfile(values, p, "ccp");
    t2 = clock();    
    std::cout << "ccp: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }
  
  return;
}

void getMpbl(const std::string mpdir, std::vector<SeqStats> &chr)
{
  std::string lineStr;
  std::vector<std::string> v;
  std::string mpfile = mpdir + "/map_fragL150_genome.txt";
  std::ifstream in(mpfile);
  if(!in) PRINTERR("Could nome open " << mpfile << ".");
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    for(auto &x: chr) {
      if(x.name == rmchr(v[0])) x.len_mpbl = stoi(v[1]);
    }
  }
  return;
}

void getMpbltable(const std::string mptable, std::vector<SeqStats> &chr)
{
  std::string lineStr;
  std::vector<std::string> v;
  std::ifstream in(mptable);
  if(!in) PRINTERR("Could nome open " << mptable << ".");
  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    for(auto &x: chr) {
      if(x.name == rmchr(v[0])) x.len_mpbl = stoi(v[1]);
    }
  }
  return;
}
