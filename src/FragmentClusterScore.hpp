/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _FRAGMENTCLUSTERSCORE_H_
#define _FRAGMENTCLUSTERSCORE_H_

#include <fstream>
#include <boost/bind.hpp>
#include "../common/inline.hpp"
#include "../common/BoostOptions.hpp"
#include "../common/seq.hpp"

namespace {
  const int32_t sizeOfvNeighborFrag(5000);
  const std::vector<int32_t> v4pnf{50, 100, 150, 500, 1000, 2000, 3000, 10000, 100000, 1000000};
}

class SeqStats;
class SeqStatsGenome;

std::vector<int8_t> genVector4FixedReadsNum(const SeqStats &chr, const double r4cmp, int32_t &numUsed4FCS, Strand::Strand strand);

class FCSstats {
  MyOpt::Opts opt;

  int32_t num4fcs;
  int32_t ng_from_fcs, ng_to_fcs, ng_step_fcs;
  double fcsread;
  double fcsflen;
  double fcs1k;
  double fcs10k;
  double fcs100k;

public:
  FCSstats():
    opt("Fragment cluster score",100),
    fcsread(0), fcsflen(0), fcs1k(0), fcs10k(0), fcs100k(0)
  {
    opt.add_options()
      ("ng_from_fcs",
       boost::program_options::value<int32_t>()->default_value(NUM_100K)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_from_fcs")),
       "fcs start of background")
      ("ng_to_fcs",
       boost::program_options::value<int32_t>()->default_value(NUM_1M)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_to_fcs")),
       "fcs end of background")
      ("ng_step_fcs",
       boost::program_options::value<int32_t>()->default_value(NUM_100K)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_step_fcs")),
       "fcs step on of background")
      ;
  }
  
  void setOpts(MyOpt::Opts &allopts) {
    allopts.add(opt);
  }
  
  void setValues(const MyOpt::Variables &values) {
    DEBUGprint("FCSstats setValues...");
    
    num4fcs     = MyOpt::getVal<int32_t>(values, "num4ssp");
    ng_from_fcs = MyOpt::getVal<int32_t>(values, "ng_from_fcs");
    ng_to_fcs   = MyOpt::getVal<int32_t>(values, "ng_to_fcs");
    ng_step_fcs = MyOpt::getVal<int32_t>(values, "ng_step_fcs");
    
    DEBUGprint("FCSstats setValues done.");
  }
    void dump()
  {
    std::cout << boost::format("FCS background region: [%d,%d], step %d\n") % ng_from_fcs % ng_to_fcs % ng_step_fcs;
    std::cout << boost::format("Read number for FCS: %d\n") % num4fcs;
  }
  int32_t getnum4fcs()   const { return num4fcs; }
  int32_t getNgFromFCS() const { return ng_from_fcs; }
  int32_t getNgToFCS()   const { return ng_to_fcs; }
  int32_t getNgStepFCS() const { return ng_step_fcs; }
  void setfcsread(const double c) { fcsread = c; }
  void setfcsflen(const double c) { fcsflen = c; }
  void setfcs1k(const double c)   { fcs1k   = c; }
  void setfcs10k(const double c)  { fcs10k  = c; }
  void setfcs100k(const double c) { fcs100k = c; }

  void printhead(std::ofstream &out) {
    out << "FCS(read)\tFCS(flen)\tFCS(1k)\tFCS(10k)\tFCS(100k)" << std::endl;
  }
  void print(std::ofstream &out) {
    out << fcsread << "\t" << fcsflen << "\t" << fcs1k << "\t"
	<< fcs10k  << "\t" << fcs100k << std::endl;
  }
};


class PropNeighborFrag {
  double sumOfvNeighborFrag;
  std::vector<int32_t> vNeighborFrag;

 public:
 PropNeighborFrag():
   sumOfvNeighborFrag(0),
   vNeighborFrag(sizeOfvNeighborFrag, 0) {}

  void setNeighborFrag(const int32_t flen, const int32_t end, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev);

  double getPNF(const int32_t i) const {
    if(i<0 || i>= sizeOfvNeighborFrag) {
      std::cerr << "error: invalid num " << i << "for getPNF max: " << sizeOfvNeighborFrag << std::endl;
      return -1;
    }
    else {
      return getratio(vNeighborFrag[i], sumOfvNeighborFrag);
    }
  }
  double getCumulativePNF(const int32_t i) const {
    if(i<0 || i>= sizeOfvNeighborFrag) {
      std::cerr << "error: invalid num " << i << "for getCumulativePNF max: " << sizeOfvNeighborFrag << std::endl;
      return -1;
    }
    else {
      double cpnf(0);
      for(int32_t j=0; j<=i; ++j) cpnf += getPNF(j);
      return cpnf;
    }
  }
};


class shiftFragVar {
  std::map<int32_t, double> mpFCS;
  std::map<int32_t, PropNeighborFrag> pnf;
  std::map<int32_t, PropNeighborFrag> pnfbg;
  int32_t lenF3;
  int32_t flen;
  double r4cmp;
  int32_t numUsed4FCS;
  bool lackOfReads;
  int32_t ng_from_fcs;
  int32_t ng_to_fcs;
  int32_t ng_step_fcs;

  long nread;
  
 public:
  shiftFragVar(const FCSstats &fcsst, const SeqStatsGenome &genome);

  void execchr(const SeqStats &chr);
  uint32_t getnumUsed4FCS() const { return numUsed4FCS; }

  double getpnfbg(const int32_t i) const {
    double v(0);
    for(auto itr = pnfbg.begin(); itr != pnfbg.end(); ++itr) {
      v += itr->second.getCumulativePNF(i);
    }
    return v/pnfbg.size();
  }
  double getFCS(const int32_t len) const {
    double diffMax(0);
    for(size_t k=0; k<sizeOfvNeighborFrag; ++k) {
      diffMax = std::max(diffMax, pnf.at(len).getCumulativePNF(k) - getpnfbg(k));
    }
    return diffMax;
  }

  void calcFCS() {
    for(auto itr = pnf.begin(); itr != pnf.end(); ++itr) {
      mpFCS[itr->first] = getFCS(itr->first);
    }
  }
  
  void outputPnf(const std::string &filename) const {
    std::ofstream out(filename);

    for(auto &x: v4pnf) out << "\tPNF len"  << x;
    for(auto &x: v4pnf) out << "\tCPNF len" << x;
    out << std::endl;

    for(size_t k=0; k<sizeOfvNeighborFrag-1; ++k) {
      out << k << "\t";
      for(auto &x: v4pnf) out << pnf.at(x).getPNF(k)           << "\t";
      for(auto &x: v4pnf) out << pnf.at(x).getCumulativePNF(k) << "\t";
      out << std::endl;
    }
  }

  void setFCSstats(FCSstats &p) {
    p.setfcsread(mpFCS.at(lenF3));
    p.setfcsflen(mpFCS.at(flen));
    p.setfcs1k(mpFCS.at(1000));
    p.setfcs10k(mpFCS.at(10000));
    p.setfcs100k(mpFCS.at(100000));
  }

  void outputFCS(const std::string &filename) const {
    if(!nread) std::cerr << filename << ": no read" << std::endl;

    std::ofstream out(filename);
    std::string str("");
    if(lackOfReads) str = " (read number is insufficient)";
    out << "Read length"     << str << "\t" << mpFCS.at(lenF3) << std::endl;
    out << "Fragment length" << str << "\t" << mpFCS.at(flen)  << std::endl;
    out << "Broad (1 kbp)"   << str << "\t" << mpFCS.at(1000)  << std::endl;
    out << "Broad (10 kbp)"  << str << "\t" << mpFCS.at(10000) << std::endl;
    out << "Strand shift\tFragment cluster score" << std::endl;
    for(auto itr = mpFCS.begin(); itr != mpFCS.end(); ++itr) {
      out << itr->first << "\t" << mpFCS.at(itr->first) << std::endl;
    }
  }
};

void makeFCSProfile(FCSstats &, const SeqStatsGenome &genome, const std::string &head, const std::string &typestr);

#endif /* _FRAGMENTCLUSTERSCORE_H_ */
