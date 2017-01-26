/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _MAPFILECLASS_H_
#define _MAPFILECLASS_H_

#include <boost/thread.hpp>
#include "seq.h"
#include "macro.h"
#include "statistics.h"
#include "bpstatus.h"


#include <fstream>
class WigArray {
  std::vector<int32_t> array;
  double geta;

  WigArray(){}
  
 public:
 WigArray(const size_t num, const int32_t val): array(num, val), geta(1000.0) {}
  ~WigArray(){}

  size_t size() const { return array.size(); }
  double getval(const size_t i) const { return array[i]/geta; }
  void setval(const size_t i, const double val) { array[i] = val*geta; }
  void addval(const size_t i, const double val) { array[i] += val*geta; }
  void multipleval(const size_t i, const double val) { array[i] *= val; }

  double getPercentile(double per) const {
    return MyStatistics::getPercentile(array, per)/geta;
  }

  void outputAsWig(std::ofstream &out, const int32_t binsize) const {
    for(size_t i=0; i<array.size(); ++i) {
      if(array[i]) out << boost::format("%1%\t%2%\n") % (i*binsize+1) % (array[i]/geta);
    }
  }
  void outputAsBedGraph(std::ofstream &out, const int32_t binsize, const std::string &name, const uint64_t chrend) const {
    uint64_t e;
    for(size_t i=0; i<array.size(); ++i) {
      if(i==array.size() -1) e = chrend; else e = (i+1) * binsize;
      if(array[i]) out << boost::format("%1% %2% %3% %4%\n") % name % (i*binsize) % e % (array[i]/geta);
    }
  }
  void outputAsBinary(std::ofstream &out) const {
    for(size_t i=0; i<array.size(); ++i) out.write((char *)&array[i], sizeof(int32_t));
  }
  void readBinary(std::ifstream &in, const int32_t nbin) const {
    for(int32_t i=0; i<nbin; ++i) in.read((char *)&array[i], sizeof(int32_t));
  }
};

class WigStats {
  enum{n_mpDist=20, n_wigDist=200};
  uint64_t sum;
 public:
  double ave, var, nb_p, nb_n, nb_p0; // pois_p0
  std::vector<uint64_t> mpDist;
  std::vector<uint64_t> wigDist;
  std::vector<double> pwigDist;

 WigStats(): sum(0), ave(0), var(0), nb_p(0), nb_n(0), nb_p0(0),
    mpDist(n_mpDist,0), wigDist(n_wigDist,0), pwigDist(n_wigDist,0) {}

  double getPoisson(const int32_t i) const {
    if(ave) return _getPoisson(i, ave);
    else return 0;
  }
  double getNegativeBinomial(const int32_t i) const {
    return _getNegativeBinomial(i, nb_p, nb_n);
  }
  double getZINB(const int32_t i) const {
    if(ave) return _getZINB(i, nb_p, nb_n, nb_p0);
    else return 0;
  }
  int32_t getwigDistthre() const {
    int32_t thre(9);
    uint64_t num;
    do{
      ++thre;
      num=0;
      for(int32_t i=0; i<thre; ++i) num += wigDist[i];
    } while(num < sum*0.8 && thre <n_wigDist-1);
#ifdef DEBUG
    std::cout << boost::format("\nthre %1%  (%2% / %3%)\n") % thre % num % sum;
#endif
    return thre;
  }
  void estimateParam() {
    int32_t thre = getwigDistthre();
    double par[thre+1];
    par[0] = thre;
    for(int32_t i=0; i<thre; ++i) par[i+1] = pwigDist[i];
    iterateZINB(&par, nb_p, nb_n, nb_p, nb_n, nb_p0);
  }
  void setpWigDist() {
    for(size_t i=0; i<wigDist.size(); ++i) pwigDist[i] = getratio(wigDist[i], sum);
  }
  void getWigStats(const WigArray &wigarray) {
    double num95 = wigarray.getPercentile(0.95);
    
    int32_t size = wigDist.size();
    std::vector<int32_t> ar;
    for(size_t i=0; i<wigarray.size(); ++i) {
      int32_t v(wigarray.getval(i));
      if(v<0) std::cout << sum << "xxx" << v << std::endl;
      ++sum;
      if(v < size) ++wigDist[v];
      if(v >= num95) continue;
      ar.push_back(v);
    }
    setpWigDist();

    MyStatistics::moment<int32_t> mm(ar, 0);
    ave = mm.getmean();
    var = mm.getvar();
    nb_p = var ? ave/var : 0;
    if(nb_p>=1) nb_p = 0.9;
    if(nb_p<=0) nb_p = 0.1; 
    nb_n = ave * nb_p /(1 - nb_p);

    //    std::cout << ave << "\t" << var << "\t" << nb_p << "\t" << nb_n<< std::endl;
    if(ave) estimateParam();
  }
  /*  void getWigStatsold(const std::vector<int32_t> &wigarray) {
    int32_t num95 = getPercentile(wigarray, 0.95);
    
    int32_t size = wigDist.size();
    std::vector<int32_t> ar;
    for(auto x: wigarray) {
      if( x<0) std::cout << sum << "xxx" << x << std::endl;
      ++sum;
      int32_t v(wigarray2value(x));
      if(v < size) ++wigDist[v];
      if(x >= num95) continue;
      ar.push_back(v);
    }
    setpWigDist();

    MyStatistics::moment<int32_t> mm(ar, 0);
    ave = mm.getmean();
    var = mm.getvar();
    nb_p = var ? ave/var : 0;
    if(nb_p>=1) nb_p = 0.9;
    if(nb_p<=0) nb_p = 0.1; 
    nb_n = ave * nb_p /(1 - nb_p);

    //    std::cout << ave << "\t" << var << "\t" << nb_p << "\t" << nb_n<< std::endl;
    if(ave) estimateParam();
    }*/
  void printwigDist(std::ofstream &out, const int32_t i) const {
    out << boost::format("%1%\t%2%\t") % wigDist[i] % pwigDist[i];
  }
  void addmpDist(const double p) {
    if(!my_range(p,0,1)) std::cout << "Warning: mappability " << p << " should be [0,1]" << std::endl;
    else ++mpDist[(int32_t)(p*n_mpDist)];
  }
  void addWigDist(const WigStats &x) {
    for(uint32_t i=0; i<wigDist.size(); ++i) wigDist[i] += x.wigDist[i];
    sum += x.sum;
    setpWigDist();
  }
  void printmpDist() const {
    uint64_t num = accumulate(mpDist.begin(), mpDist.end(), 0);
    for(size_t i=0; i<mpDist.size(); ++i)
      std::cout << boost::format("~%1%%%\t%2%\t%3%\n") % ((i+1)*100/mpDist.size()) % mpDist[i] % getratio(mpDist[i], num); 
  }
  void printPoispar(std::ofstream &out) const {
    out << boost::format("%1$.3f\t%2$.3f\t") % ave % var;
  }
  void printZINBpar(std::ofstream &out) const {
    out << boost::format("%1%\t%2%\t%3%") % nb_p % nb_n % nb_p0;
  }
};

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
void calcdepth(T &obj, const int32_t flen)
{
  uint64_t lenmpbl = obj.getlenmpbl();
  double d = lenmpbl ? getratio(obj.getnread_nonred(Strand::BOTH) * flen, lenmpbl): 0;
  obj.setdepth(d);
}

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
  WigStats ws;
    
 SeqStats(std::string s, int32_t l=0, int32_t binsize=0):
  name(rmchr(s)), len(l), len_mpbl(l), 
    Greekchr(false), depth(0), nread_inbed(0),
    nbp(0), ncov(0), ncovnorm(0), sizefactor(0)
    {
    nbin = binsize ? (l/binsize +1) : 0;
  }
  
  virtual ~SeqStats(){}

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
  strandData & getStrandref (const Strand::Strand strand) {
    return seq[strand];
  }
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
    for(auto x: array) {
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

  friend void getMpbl(const std::string, std::vector<SeqStats> &chr);
  friend void getMpbltable(const std::string, std::vector<SeqStats> &chr);
};


#endif /* _MAPFILECLASS_H_ */
