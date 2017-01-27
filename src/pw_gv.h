/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <fstream>
#include <numeric>
#include <boost/thread.hpp>
#include "readdata.h"
#include "mthread.h"
#include "LibraryComplexity.hpp"
#include "Mapfile.hpp"

class SSPstats {
  double nsc;
  double rsc;
  double backgroundUniformity;
  double fcsread;
  double fcsflen;
  double fcs1k;
  double fcs10k;
  double fcs100k;

 public:
 SSPstats(): nsc(0), rsc(0), backgroundUniformity(0), fcsflen(0), fcs1k(0), fcs10k(0), fcs100k(0) {}

  void setnsc(const double c) { nsc = c; }
  void setrsc(const double c) { rsc = c; }
  void setbu(const double c) { backgroundUniformity = c; }
  void setfcsread(const double c) { fcsread = c; }
  void setfcsflen(const double c) { fcsflen = c; }
  void setfcs1k(const double c) { fcs1k = c; }
  void setfcs10k(const double c) { fcs10k = c; }
  void setfcs100k(const double c) { fcs100k = c; }

  void printhead(std::ofstream &out) {
    out << "NSC\tRSC\tbackground uniformity\t"
	<< "FCS(read)\tFCS(flen)\tFCS(1k)\tFCS(10k)\tFCS(100k)"
	<< std::endl;
  }
  void print(std::ofstream &out) {
    out << nsc << "\t" << rsc << "\t" << backgroundUniformity << "\t"
	<< fcsread << "\t" << fcsflen << "\t" << fcs1k << "\t" << fcs10k
	<< "\t" << fcs100k << std::endl;
  }
};

class Mapfile: private Uncopyable {
  bool Greekchr;

  std::string samplename;
  std::string oprefix;
  std::string obinprefix;

  bool lackOfRead4GenomeCov;
  bool lackOfRead4FragmentVar;
  std::vector<Peak> vPeak;

  // GC bias
  int32_t maxGC;

  void setlchr() {
    uint64_t lenmax(0);
    for(auto itr = genome.chr.begin(); itr != genome.chr.end(); ++itr) {
      if(lenmax < itr->getlenmpbl()) {
	lenmax = itr->getlenmpbl();
	lchr = itr;
      }
    }
  }
  
 public:
  SeqStatsGenome genome;
  WigStats wsGenome;
  std::vector<SeqWigStats>::iterator lchr; // longest chromosome

  // for SSP
  SSPstats sspst;

  class LibComp complexity;
  // Wigdist
  int32_t nwigdist;
  std::vector<int32_t> wigDist;
  
 Mapfile(const MyOpt::Variables &values):
  Greekchr(false),
    samplename(values["output"].as<std::string>()),
    lackOfRead4GenomeCov(false), lackOfRead4FragmentVar(false),
    maxGC(0), genome(values), complexity(values) {
      setlchr();

      oprefix = values["odir"].as<std::string>() + "/" + values["output"].as<std::string>();
      obinprefix = oprefix + "." + IntToString(values["binsize"].as<int32_t>());
    }

  void outputSSPstats() {
    std::string filename = getprefix() + ".stats.txt";
    std::ofstream out(filename);
    out << "Sample\ttotal read number\tnonredundant read number\t"
	<< "read length\tfragment length\t";
    sspst.printhead(out);
    out << samplename << "\t" << genome.getnread(Strand::BOTH) << "\t" << genome.getnread_nonred(Strand::BOTH) << "\t"
	<< genome.dflen.getlenF3() << "\t" << genome.dflen.getflen() << "\t";
    sspst.print(out);
  }
  
  void setmaxGC(const int32_t m) { maxGC = m; }
  int32_t getmaxGC() const {return maxGC; }

  void lackOfRead4GenomeCov_on() { lackOfRead4GenomeCov = true; }
  bool islackOfRead4GenomeCov() const { return lackOfRead4GenomeCov; };
  void printPeak(const MyOpt::Variables &values) const {
    std::string filename = getbinprefix() + ".peak.xls";
    std::ofstream out(filename);

    vPeak[0].printHead(out);
    for(uint32_t i=0; i<vPeak.size(); ++i) {
      vPeak[i].print(out, i, values["binsize"].as<int32_t>());
    }
  }
  std::string getprefix() const { return oprefix; }
  std::string getbinprefix() const { return obinprefix; }

  void addPeak(const Peak &peak) {
    vPeak.push_back(peak);
  }
  void renewPeak(const int32_t i, const double val, const double p) {
    vPeak[vPeak.size()-1].renew(i, val, p);
  }

  void estimateZINB() {
    int32_t thre = wsGenome.getwigDistthre();
    double par[thre+1];
    par[0] = thre;
    for(int32_t i=0; i<thre; ++i) par[i+1] = wsGenome.pwigDist[i];

    iterateZINB(&par, lchr->ws.nb_p, lchr->ws.nb_n, wsGenome.nb_p, wsGenome.nb_n, wsGenome.nb_p0);

    return;
  }
};

#endif /* _PW_GV_H_ */
