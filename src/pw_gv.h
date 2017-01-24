/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <fstream>
#include <numeric>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include "macro.h"
#include "util.h"
#include "readdata.h"
#include "statistics.h"

namespace MyOpt {
  using Variables = boost::program_options::variables_map;
  using Opts      = boost::program_options::options_description;
}

void printDist(std::ofstream &out, const std::vector<int32_t> v, const std::string str, const uint64_t nread);

class Fragment {
public:
  std::string chr;
  int32_t F3;
  Strand strand;
  int32_t fraglen;
  int32_t readlen_F3;

 Fragment(): fraglen(0), readlen_F3(0) {}
 void addSAM(const std::vector<std::string> &v, const bool pair, const int32_t sv) {
   chr = rmchr(v[2]);
   readlen_F3 = v[9].length();
   if(pair) fraglen = abs(stoi(v[8]));
   else fraglen = readlen_F3;
   if(sv&16) {
     strand = STRAND_MINUS;
     F3 = stoi(v[3]) + readlen_F3 -1;
   } else {
     strand = STRAND_PLUS;
     F3 = stoi(v[3]) -1;
   }
 }
 void print() const {
   BPRINT("chr:%1%\tposi:%2%\tstrand:%3%\tfraglen:%4%\treadlen:%5%\n") % chr % F3 % strand % fraglen % readlen_F3;
  }
};

class Read {
  int32_t weight;
  enum {WeightNum=1000};
 public:
  int32_t F3;
  int32_t F5;
  int32_t duplicate;
  int32_t inpeak;
  
 Read(const Fragment &frag): weight(WeightNum), F3(frag.F3), duplicate(0), inpeak(0) {
    if(frag.strand == STRAND_PLUS) F5 = frag.F3 + frag.fraglen;
    else F5 = frag.F3 - frag.fraglen;
  }
  double getWeight() const {
    return weight/static_cast<double>(WeightNum);
  }
  void multiplyWeight(const double w) { weight *= w; }
};

class strandData {
 public:
  std::vector<Read> vRead;
  uint64_t nread;
  uint64_t nread_nonred;
  uint64_t nread_red;
  double nread_rpm;
  double nread_afterGC;

 strandData(): nread(0), nread_nonred(0), nread_red(0), nread_rpm(0), nread_afterGC(0) {}
  void setnread() { nread = vRead.size(); }
  void print() const {
    std::cout << nread << "\t" << nread_nonred << "\t" << nread_red << "\t" << nread_rpm << "\t" << nread_afterGC << std::endl;
  }
  void printnonred(std::ofstream &out)  const { printr(out, nread_nonred,  nread); }
  void printred(std::ofstream &out)     const { printr(out, nread_red,     nread); }
  void printafterGC(std::ofstream &out) const { printr(out, nread_afterGC, nread); }
  void addReadAfterGC(const double w, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    nread_afterGC += w;
  }
  void setnread_nonread_nofilter() {
    nread_nonred = nread;
  }
};

class Wigstats {
  enum{n_mpDist=20, n_wigDist=200};
  uint64_t sum;
 public:
  double ave, var, nb_p, nb_n, nb_p0; // pois_p0
  std::vector<uint64_t> mpDist;
  std::vector<uint64_t> wigDist;
  std::vector<double> pwigDist;

 Wigstats(): sum(0), ave(0), var(0), nb_p(0), nb_n(0), nb_p0(0),
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
    BPRINT("\nthre %1%  (%2% / %3%)\n") % thre % num % sum;
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
    for(size_t i=0; i<wigDist.size(); ++i) pwigDist[i] = wigDist[i]/static_cast<double>(sum);
  }
  void getWigStats(const std::vector<int32_t> &wigarray) {
    int32_t num95 = getPercentile(wigarray, 0.95);
    
    int32_t size = wigDist.size();
    std::vector<int32_t> ar;
    for(auto x: wigarray) {
    if( x<0) std::cout << sum << "xxx" << x << std::endl;
    ++sum;
      int32_t v = WIGARRAY2VALUE(x);
      if(v < size) ++wigDist[v];
      if(x >= num95) continue;
      ar.push_back(v);
    }
    setpWigDist();

    moment<int32_t> mm(ar, 0);
    ave = mm.getmean();
    var = mm.getvar();
    nb_p = var ? ave/var : 0;
    if(nb_p>=1) nb_p = 0.9;
    if(nb_p<=0) nb_p = 0.1; 
    nb_n = ave * nb_p /(1 - nb_p);

    //    std::cout << ave << "\t" << var << "\t" << nb_p << "\t" << nb_n<< std::endl;
    if(ave) estimateParam();
  }
  void printwigDist(std::ofstream &out, const int32_t i) const {
    out << boost::format("%1%\t%2%\t") % wigDist[i] % pwigDist[i];
  }
  void addmpDist(const double p) {
    if(!my_range(p,0,1)) std::cout << "Warning: mappability " << p << " should be [0,1]" << std::endl;
    else ++mpDist[(int32_t)(p*n_mpDist)];
  }
  void addWigDist(const Wigstats &x) {
    for(uint32_t i=0; i<wigDist.size(); ++i) wigDist[i] += x.wigDist[i];
    sum += x.sum;
    setpWigDist();
  }
  void printmpDist() const {
    uint64_t num = accumulate(mpDist.begin(), mpDist.end(), 0);
    for(size_t i=0; i<mpDist.size(); ++i)
      BPRINT("~%1%%%\t%2%\t%3%\n") % ((i+1)*100/mpDist.size()) % mpDist[i] % (mpDist[i]/static_cast<double>(num)); 
  }
  void printPoispar(std::ofstream &out) const {
    out << boost::format("%1$.3f\t%2$.3f\t") % ave % var;
  }
  void printZINBpar(std::ofstream &out) const {
    out << boost::format("%1%\t%2%\t%3%") % nb_p % nb_n % nb_p0;
  }
};

class SeqStats {
 protected:
  bool yeast;
  uint64_t len, len_mpbl;
  double weight4rpm;
  /* FRiP */
  uint64_t nread_inbed;
  //  GenomeCoverage gcov;
  uint64_t nbp, ncov, ncovnorm;
  double gcovRaw, gcovNorm;
 public:
  std::string name;
  int32_t nbin;

  Wigstats ws;
  
  strandData seq[STRANDNUM];
  double depth;
  
 SeqStats(std::string s, int32_t l=0): yeast(false), len(l), len_mpbl(l), weight4rpm(0), nread_inbed(0), nbp(0), ncov(0), ncovnorm(0), gcovRaw(0), gcovNorm(0) , nbin(0), depth(0) {
    name = rmchr(s);
  }
  virtual ~SeqStats(){}
  uint64_t getNreadInbed() const { return nread_inbed; }
  void addfrag(const Fragment &frag) {
    Read r(frag);
    seq[frag.strand].vRead.push_back(r);
  }
  uint64_t bothnread ()         const { return seq[STRAND_PLUS].nread         + seq[STRAND_MINUS].nread; }
  uint64_t bothnread_nonred ()  const { return seq[STRAND_PLUS].nread_nonred  + seq[STRAND_MINUS].nread_nonred; }
  uint64_t bothnread_red ()     const { return seq[STRAND_PLUS].nread_red     + seq[STRAND_MINUS].nread_red; }
  uint64_t bothnread_rpm ()     const { return seq[STRAND_PLUS].nread_rpm     + seq[STRAND_MINUS].nread_rpm; }
  uint64_t bothnread_afterGC () const { return seq[STRAND_PLUS].nread_afterGC + seq[STRAND_MINUS].nread_afterGC; }

  uint64_t getlen()      const { return len; }
  uint64_t getlenmpbl()  const { return len_mpbl; }
  uint64_t getnbp()      const { return nbp; }
  uint64_t getncov()     const { return ncov; }
  uint64_t getncovnorm() const { return ncovnorm; }
  double getpmpbl()  const { return static_cast<double>(len_mpbl)/len; }
  int32_t getnbin()      const { return nbin; }
  double getweight4rpm() const { return weight4rpm; }
  void print() const {
    std::cout << name << "\t" << len << "\t" << len_mpbl << "\t" << bothnread() << "\t" << bothnread_nonred() << "\t" << bothnread_red() << "\t" << bothnread_rpm() << "\t" << bothnread_afterGC()<< "\t" << depth << std::endl;
  }
  void printGcov(std::ofstream &out, const bool lackOfRead4GenomeCov) const {
      if(lackOfRead4GenomeCov) out << boost::format("%1$.3f\t(%2$.3f)\t") % gcovRaw % gcovNorm;
      else out << boost::format("%1$.3f\t%2$.3f\t")   % gcovRaw % gcovNorm;
  }
  void calcdepth(const int32_t flen) {
    depth = len_mpbl ? bothnread_nonred() * flen / static_cast<double>(len_mpbl): 0;
  }
  void setF5(int32_t flen) {
    int32_t d;
    for(int32_t strand=0; strand<STRANDNUM; ++strand) {
      if(strand == STRAND_PLUS) d = flen; else d = -flen;
      for(auto &x: seq[strand].vRead) x.F5 = x.F3 + d;
    }
  }
  double getFRiP() const {
    return nread_inbed/static_cast<double>(bothnread_nonred());
  }
  void setWeight(const double weight) {
    weight4rpm = weight;
    for(int32_t i=0; i<STRANDNUM; i++) seq[i].nread_rpm = seq[i].nread_nonred * weight4rpm;
  }
  void calcGcov(const std::vector<int8_t> &array) {
    for(auto x:array) {
      if(x >= MAPPABLE)     ++nbp;      // MAPPABLE || COVREAD_ALL || COVREAD_NORM
      if(x >= COVREAD_ALL)  ++ncov;     // COVREAD_ALL || COVREAD_NORM
      if(x == COVREAD_NORM) ++ncovnorm;
    }
    gcovRaw  = nbp ? ncov     / static_cast<double>(nbp): 0;
    gcovNorm = nbp ? ncovnorm / static_cast<double>(nbp): 0;
  }
  void yeaston() { yeast = true; }

  bool isautosome() const {
    int32_t chrnum(0);
    try {
      chrnum = stoi(name);
    } catch (std::invalid_argument e) {  // 数値以外
      if(yeast) { 
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
  friend void calcFRiP(SeqStats &, const std::vector<bed>);
};

class SeqStatsGenome: public SeqStats {
  std::vector<bed> vbed;
  void readGenomeTable(const std::string &, const int32_t);
  
  void printstats() const {
    std::cout << "name\tlength\tlen_mpbl\tread num\tnonred num\tred num\tnormed\tafterGC\tdepth" << std::endl;
    print();
    for(auto x:chr) x.print();
  }
  
  std::vector<sepchr> getVsepchr(const int32_t numthreads){
    std::vector<sepchr> vsep;
    
    uint32_t sepsize = len/numthreads;
    for(uint32_t i=0; i<chr.size(); ++i) {
      uint32_t s = i;
      uint64_t len(0);
      while(len < sepsize && i<chr.size()) {
	len += chr[i].getlen();
	i++;
      }
      i--;
      uint32_t e = i;
      sepchr sep(s,e);
      vsep.push_back(sep);
    }
    return vsep;
  }
  
 public:
  std::vector<SeqStats> chr;
  std::vector<sepchr> vsepchr;
  
 SeqStatsGenome(const MyOpt::Variables &values): SeqStats("Genome") {
    
    readGenomeTable(values["gt"].as<std::string>(), values["binsize"].as<int32_t>());
    if(values.count("mp")) getMpbl(values["mp"].as<std::string>(), chr);
    else if(values.count("mptable")) getMpbltable(values["mptable"].as<std::string>(), chr);
    for(auto &x:chr) {
      len      += x.getlen();
      //     len_mpbl += x.getlenmpbl();
      // nbin     += x.getnbin();
    }
  
    // yeast
    for(auto x:chr) if(x.name == "I" || x.name == "II" || x.name == "III") yeast = true;
    for(auto &x:chr) if(yeast) x.yeaston();
    
    // sepchr
    vsepchr = getVsepchr(values["threads"].as<int32_t>());
    
#ifdef DEBUG
    std::cout << "chr\tautosome" << std::endl;
    for(auto x:chr) std::cout << x.name << "\t" << x.isautosome() << std::endl;
    for(uint32_t i=0; i<vsepchr.size(); i++)
      std::cout << "thread " << (i+1) << ": "
		<< vsepchr[i].s << "-" << vsepchr[i].e
		<< std::endl;
    printstats();
#endif
  }

  uint64_t getlenmpbl() const {
    uint64_t len_mpbl(0);
    for(auto &x:chr) len_mpbl += x.getlenmpbl();
    return len_mpbl;
  }
  int32_t getnbin() const {
    int32_t nbin(0);
    for(auto &x:chr) nbin += x.getnbin();
    return nbin;
  }
  double getpmpbl() const { return static_cast<double>(getlenmpbl())/len; }

  void setnread() {
    for(auto &x:chr) {
      for(int32_t i=0; i<STRANDNUM; i++) {
	x.seq[i].setnread();
	seq[i].nread += x.seq[i].nread;
      }
    }
  }
  void setnread_red() {
    for(auto &x:chr) {
      for(int32_t i=0; i<STRANDNUM; i++) {
	seq[i].nread_nonred += x.seq[i].nread_nonred;
	seq[i].nread_red    += x.seq[i].nread_red;
      }
    }
  }
  void setnread2nread_red() {
    for(auto &x:chr) {
      for(int32_t i=0; i<STRANDNUM; i++) x.seq[i].setnread_nonread_nofilter();
    }
  }
  void setbed(const std::string bedfilename) {
    isFile(bedfilename);
    vbed = parseBed<bed>(bedfilename);
    //    printBed(vbed);
  }
  void setFRiP() {
    std::cout << "calculate FRiP score.." << std::flush;
    for(auto &x: chr) {
      calcFRiP(x, vbed);
      nread_inbed += x.getNreadInbed();
    }
    std::cout << "done." << std::endl;
    return;
  }
   void calcdepthAll(const int32_t flen) {
     for (auto &x:chr) x.calcdepth(flen);
     calcdepth(flen);
  }
   void setF5All(const int32_t flen) {
    for (auto &x:chr) x.setF5(flen);
    setF5(flen);
  }
  std::vector<bed> getvbed() const { return vbed; }
  void addGcov(const int32_t i, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    nbp      += chr[i].getnbp();
    ncov     += chr[i].getncov();
    ncovnorm += chr[i].getncovnorm();
    gcovRaw  = nbp ? ncov / static_cast<double>(nbp): 0;
    gcovNorm = nbp ? ncovnorm / static_cast<double>(nbp): 0;
  }
};

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
  bool yeast;
  const int32_t ReadMax=200;
  const int32_t FragMax=1000;
  int32_t lenF3;
  int32_t lenF5;
  int32_t eflen;
  int32_t flen_def;
  std::vector<int32_t> vlenF3;
  std::vector<int32_t> vlenF5;
  std::vector<int32_t> vflen;

  std::string samplename;
  std::string oprefix;
  std::string obinprefix;
  
  // PCR bias
  int32_t thre4filtering;
  int32_t nt_all, nt_nonred, nt_red;
  bool lackOfRead4Complexity;
  bool lackOfRead4GenomeCov;
  bool lackOfRead4FragmentVar;
  double r4cmp;
  std::vector<Peak> vPeak;

  // GC bias
  int32_t maxGC;

  // for SSP
  SSPstats sspst;
  
 public:
  SeqStatsGenome genome;
  std::vector<SeqStats>::iterator lchr; // longest chromosome

  // Wigdist
  int32_t nwigdist;
  std::vector<int32_t> wigDist;
  
 Mapfile(const MyOpt::Variables &values):
  yeast(false),
    lenF3(0), lenF5(0), eflen(0), flen_def(values["flen"].as<int32_t>()),
    vlenF3(ReadMax,0), vlenF5(ReadMax,0), vflen(FragMax,0),
    samplename(values["output"].as<std::string>()),
    thre4filtering(0), nt_all(0), nt_nonred(0), nt_red(0),
    lackOfRead4Complexity(false), lackOfRead4GenomeCov(false), lackOfRead4FragmentVar(false), r4cmp(0), maxGC(0), genome(values)
    {
      uint64_t lenmax(0);
      for(auto itr = genome.chr.begin(); itr != genome.chr.end(); ++itr) {
	if(lenmax < itr->getlenmpbl()) {
	  lenmax = itr->getlenmpbl();
	  lchr = itr;
	}
      }

      oprefix = values["odir"].as<std::string>() + "/" + values["output"].as<std::string>();
      obinprefix = oprefix + "." + IntToString(values["binsize"].as<int32_t>());
    }

  void setSSPstats(const double bu, const double nsc, const double rsc) {
    sspst.setnsc(nsc);
    sspst.setrsc(rsc);
    sspst.setbu(bu);
  }
  void setFCSstats(const double fcsread, const double fcsflen, const double fcs1k, const double fcs10k, const double fcs100k) {
    sspst.setfcsread(fcsread);
    sspst.setfcsflen(fcsflen);
    sspst.setfcs1k(fcs1k);
    sspst.setfcs10k(fcs10k);
    sspst.setfcs100k(fcs100k);
  }
  void printSSPstats() {
    std::string filename = getprefix() + ".stats.txt";
    std::ofstream out(filename);
    out << "Sample\ttotal read number\tnonredundant read number\t"
	<< "read length\tfragment length\t";
    sspst.printhead(out);
    out << samplename << "\t" << genome.bothnread() << "\t" << genome.bothnread_nonred() << "\t"
	<< lenF3 << "\t" << eflen << "\t";
    sspst.print(out);
  }
  
  void setmaxGC(const int32_t m) { maxGC = m; }
  int32_t getmaxGC() const {return maxGC; }
  void readGenomeTable(const MyOpt::Variables &values);

  void lackOfRead4Complexity_on() { lackOfRead4Complexity = true; }
  void lackOfRead4GenomeCov_on() { lackOfRead4GenomeCov = true; }
  bool islackOfRead4GenomeCov() const { return lackOfRead4GenomeCov; };
  void setthre4filtering(const MyOpt::Variables &values) {
    if(values["thre_pb"].as<int32_t>()) thre4filtering = values["thre_pb"].as<int32_t>();
    else thre4filtering = std::max(1, (int32_t)(genome.bothnread() *10/static_cast<double>(genome.getlenmpbl())));
    std::cout << "Checking redundant reads: redundancy threshold " << thre4filtering << std::endl;
  }
  int32_t getthre4filtering() const { return thre4filtering; };
  void setr4cmp(const double r) { r4cmp = r; }
  double getr4cmp() const { return r4cmp; }
  void incNtAll() { ++nt_all; }
  void incNtNonred() { ++nt_nonred; }
  void incNtRed() { ++nt_red; }
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
  void setFraglen(const MyOpt::Variables &values) {
    lenF3 = getmaxi(vlenF3);
    if(values.count("pair")) {
      lenF5 = getmaxi(vlenF5);
      eflen = getmaxi(vflen);
    }
  }
  void printFlen(const MyOpt::Variables &values, std::ofstream &out) const {
    if(!values.count("nomodel")) out << "Estimated fragment length: " << eflen << std::endl;
    else out << "Predefined fragment length: " << flen_def << std::endl;
  }
  void addF5(const int32_t readlen_F5) { ++vlenF5[readlen_F5]; }
  void addfrag(const Fragment &frag) {
    ++vlenF3[frag.readlen_F3];
    ++vflen[frag.fraglen];
    int8_t on(0);
    for(auto &x:genome.chr) {
      if(x.name == frag.chr) {
	x.addfrag(frag);
	on++;
      }
    }
    if(!on) std::cerr << "Warning: " << frag.chr << " is not in genometable." << std::endl;
  }

  /*  void setbed(const std::string bedfilename) {
    isFile(bedfilename);
    vbed = parseBed<bed>(bedfilename);
    //    printBed(vbed);
    }*/
  void outputDistFile(const MyOpt::Variables &values)
  {
    std::string outputfile = oprefix + ".readlength_dist.csv";
    std::ofstream out(outputfile);
    printDist(out, vlenF3, "F3", genome.bothnread());
    if(values.count("pair")) printDist(out, vlenF5, "F5", genome.bothnread());
    out.close();
    
    if(values.count("pair")) {
      outputfile = oprefix + ".fraglen_dist.xls";
      std::ofstream out(outputfile);
      printDist(out, vflen, "Fragmemt", genome.bothnread());
    }
  }

  void printComplexity(std::ofstream &out) const {
    if(lackOfRead4Complexity) out << boost::format("Library complexity: (%1$.3f) (%2%/%3%)\n") % complexity() % nt_nonred % nt_all;
    else out << boost::format("Library complexity: %1$.3f (%2%/%3%)\n") % complexity() % nt_nonred % nt_all;
  }
  double complexity() const { return nt_nonred/static_cast<double>(nt_all); }
  void printstats() const {
    std::cout << "name\tlength\tlen_mpbl\tread num\tnonred num\tred num\tnormed\tafterGC\tdepth" << std::endl;
    genome.print();
    for (auto x:genome.chr) x.print();
  }
  void seteflen(const int32_t len) { eflen = len; }
  int32_t getlenF3() const { return lenF3; }
  int32_t getflen(const MyOpt::Variables &values) const {
    int32_t flen;
    if(!values.count("nomodel") || values.count("pair")) flen = eflen;
    else flen = flen_def;
    return flen;
  }
  void addPeak(const Peak &peak) {
    vPeak.push_back(peak);
  }
  void renewPeak(const int32_t i, const double val, const double p) {
    vPeak[vPeak.size()-1].renew(i, val, p);
  }

  void estimateZINB() {
    int32_t thre = genome.ws.getwigDistthre();
    double par[thre+1];
    par[0] = thre;
    for(int32_t i=0; i<thre; ++i) par[i+1] = genome.ws.pwigDist[i];

    //    iteratePoisson(&par, lchr->ave, genome.ave, genome.pois_p0);
    iterateZINB(&par, lchr->ws.nb_p, lchr->ws.nb_n, genome.ws.nb_p, genome.ws.nb_n, genome.ws.nb_p0);

    return;
  }
};

#endif /* _PW_GV_H_ */
