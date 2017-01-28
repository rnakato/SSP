/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _SEQSTATS_HPP_
#define _SEQSTATS_HPP_

#include <boost/thread.hpp>
#include <boost/algorithm/string.hpp> 
#include "seq.h"
#include "macro.h"
#include "statistics.h"
#include "bpstatus.h"

class strandData {
 public:
  std::vector<Read> vRead;
  uint64_t nread_nonred;
  uint64_t nread_red;
  double nread_rpm;
  double nread_afterGC;

 strandData(): nread_nonred(0), nread_red(0), nread_rpm(0), nread_afterGC(0) {}
  size_t getnread() const { return vRead.size(); }
  void print() const {
    std::cout << getnread() << "\t" << nread_nonred << "\t" << nread_red << "\t" << nread_rpm << "\t" << nread_afterGC << std::endl;
  }

  void setnread_nonread_nofilter() {
    nread_nonred = getnread();
  }
};

template <class T>
void printSeqStats(const T &obj)
{
  std::cout << obj.getname() << "\t" << obj.getlen() << "\t" << obj.getlenmpbl() << "\t"
	    << obj.getnread(Strand::BOTH)        << "\t"
	    << obj.getnread_nonred(Strand::BOTH) << "\t"
	    << obj.getnread_red(Strand::BOTH)    << "\t"
	    << obj.getnread_rpm(Strand::BOTH)    << "\t"
	    << obj.getnread_afterGC(Strand::BOTH)<< "\t"
	    << obj.getdepth() << std::endl;
}

class SeqStats {
  enum {STRANDNUM=2};

  std::string name;
  uint64_t len, len_mpbl;
  int32_t nbin;
  strandData seq[STRANDNUM];
  bool Greekchr;
  double depth;
  uint64_t nread_inbed;
  uint64_t nbp, ncov, ncovnorm;
  double sizefactor;

 public:    
 SeqStats(std::string s, int32_t l=0, int32_t binsize=0):
  name(rmchr(s)), len(l), len_mpbl(l), 
  Greekchr(false), depth(0), nread_inbed(0),
  nbp(0), ncov(0), ncovnorm(0), sizefactor(0)
  {
    nbin = binsize ? (l/binsize +1) : 0;
  }

  std::string getname() const { return name; }
  void addfrag(const Fragment &frag) {
    Read r(frag);
    seq[frag.strand].vRead.push_back(r);
  }
  uint64_t getnread (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].getnread() + seq[Strand::REV].getnread();
    else return seq[strand].getnread();
  }
  uint64_t getnread_nonred (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].nread_nonred + seq[Strand::REV].nread_nonred;
    else return seq[strand].nread_nonred;
  }
  uint64_t getnread_red (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].nread_red + seq[Strand::REV].nread_red;
    else return seq[strand].nread_red;
  }
  uint64_t getnread_rpm (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].nread_rpm + seq[Strand::REV].nread_rpm;
    else return seq[strand].nread_rpm;
  }
  uint64_t getnread_afterGC (const Strand::Strand strand) const {
    if(strand==Strand::BOTH) return seq[Strand::FWD].nread_afterGC + seq[Strand::REV].nread_afterGC;
    else return seq[strand].nread_afterGC;
  }
  uint64_t getnread_inbed() const { return nread_inbed; }

  void incNReadNonred (const Strand::Strand strand) { ++seq[strand].nread_nonred; }
  void incNReadRed    (const Strand::Strand strand) { ++seq[strand].nread_red;    }

  const std::vector<Read> & getvReadref (const Strand::Strand strand) const {
    return seq[strand].vRead;
  }
  std::vector<Read> & getvReadref_notconst (const Strand::Strand strand) {
    return seq[strand].vRead;
  }
  void setdepth(const double d) { depth = d; }
  void addReadAfterGC(const Strand::Strand strand, const double w, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    seq[strand].nread_afterGC += w;
  }
  
  uint64_t getlen()       const { return len; }
  uint64_t getlenmpbl()   const { return len_mpbl; }
  double   getpmpbl()     const { return getratio(getlenmpbl(), getlen()); }
  uint64_t getnbp()       const { return nbp; }
  uint64_t getncov()      const { return ncov; }
  uint64_t getncovnorm()  const { return ncovnorm; }
  int32_t  getnbin()      const { return nbin; }
  double   getsizefactor()const { return sizefactor; }
  double   getdepth()     const { return depth; }

  void setF5(int32_t flen) {
    int32_t d;
    for (auto strand: {Strand::FWD, Strand::REV}) {
      if(strand == Strand::FWD) d = flen; else d = -flen;
      for(auto &x: seq[strand].vRead) x.F5 = x.F3 + d;
    }
  }
  void setFRiP(const uint64_t n) { nread_inbed = n; }
  double getFRiP() const {
    return getratio(nread_inbed, getnread_nonred(Strand::BOTH));
  }
  void setsizefactor(const double w) {
    sizefactor = w;
    for (auto strand: {Strand::FWD, Strand::REV}) seq[strand].nread_rpm = seq[strand].nread_nonred * sizefactor;
  }
  void calcGcov(const std::vector<BpStatus> &array) {
    for(auto &x: array) {
      if(x >= BpStatus::MAPPABLE)     ++nbp;      // MAPPABLE || COVREAD_ALL || COVREAD_NORM
      if(x >= BpStatus::COVREAD_ALL)  ++ncov;     // COVREAD_ALL || COVREAD_NORM
      if(x == BpStatus::COVREAD_NORM) ++ncovnorm;
    }
  }
  void Greekchron() { Greekchr = true; }

  bool isautosome() const {
    int32_t chrnum(0);
    try {
      chrnum = stoi(name);
    } catch (std::invalid_argument e) {  // 数値以外
      if(Greekchr) { 
	if(name=="I")         chrnum = 1;
	else if(name=="II")   chrnum = 2;
	else if(name=="III")  chrnum = 3;
	else if(name=="IV")   chrnum = 4;
	else if(name=="V")    chrnum = 5;
	else if(name=="VI")   chrnum = 6;
	else if(name=="VII")  chrnum = 7;
	else if(name=="VIII") chrnum = 8;
	else if(name=="IX")   chrnum = 9;
	else if(name=="X")    chrnum = 10;
	else if(name=="XI")   chrnum = 11;
	else if(name=="XII")  chrnum = 12;
	else if(name=="XIII") chrnum = 13;
	else if(name=="XIV")  chrnum = 14;
	else if(name=="XV")   chrnum = 15;
	else if(name=="XVI")  chrnum = 16;
      }
      if(name=="2L") chrnum = 1;
      if(name=="2R") chrnum = 2;
      if(name=="3L") chrnum = 3;
      if(name=="3R") chrnum = 4;
    }
    if(chrnum) return true;
    else       return false;
  }

  void getMptable(const std::string &mptable)
  {
    std::string lineStr;
    std::vector<std::string> v;
    std::ifstream in(mptable);
    if(!in) PRINTERR("Could nome open " << mptable << ".");
    while (!in.eof()) {
      getline(in, lineStr);
      if(lineStr.empty() || lineStr[0] == '#') continue;
      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      if(name == rmchr(v[0])) len_mpbl = stoi(v[1]);
    }
    return;
  }
};

#endif /* _MAPFILECLASS_H_ */
