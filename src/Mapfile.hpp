/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _MAPFILE_HPP_
#define _MAPFILE_HPP_

#include <boost/format.hpp>
#include "BoostOptions.hpp"
#include "SeqStats.hpp"
#include "wigstats.h"
#include "mthread.h"
#include "readdata.h"
#include "util.h"

class FragmentLengthDist {
  enum {ReadLenMax=200, FragLenMax=1000};
  int32_t flen_ssp;
  int32_t flen_def;
  std::vector<int32_t> vlenF3;
  std::vector<int32_t> vlenF5;
  std::vector<int32_t> vflen;
  bool nomodel;
  bool pairedend;

  template <class T>
  void printVector(std::ofstream &out, const std::vector<T> v, const std::string &str, const uint64_t nread)
  {
    out << str << " length distribution" << std::endl;
    out << "length\tnumber\tproportion" << std::endl;
    for(size_t i=0; i<v.size(); ++i)
      if(v[i]) out << boost::format("%1%\t%2%\t%3%\n") % i % v[i] % getratio(v[i], nread);
  }
 
 public:
 FragmentLengthDist(const MyOpt::Variables &values):
  flen_ssp(0),
    flen_def(values["flen"].as<int32_t>()),
    vlenF3(ReadLenMax,0),
    vlenF5(ReadLenMax,0),
    vflen(FragLenMax,0),
    nomodel(values.count("nomodel")),
    pairedend(values.count("pair")) {}

  int32_t getlenF3 () const { return getmaxi(vlenF3); }
  int32_t getlenF5 () const { return getmaxi(vlenF5); }
  int32_t getflen4paired () const { return getmaxi(vflen); }

  void setflen_ssp(const int32_t len) { flen_ssp = len; }
  int32_t getflen() const {
    if(pairedend)     return getflen4paired();
    else if(!nomodel) return flen_ssp;
    else              return flen_def;
  }
  void printFlen(std::ofstream &out) const {
    if(pairedend)     out << "Most likely fragment length: " << getflen4paired() << std::endl;
    else if(!nomodel) out << "Estimated fragment length: "   << flen_ssp << std::endl;
    else              out << "Predefined fragment length: "  << flen_def << std::endl;
  }
  void addF3(const int32_t lenF3) { ++vlenF3[lenF3]; }
  void addF5(const int32_t lenF5) { ++vlenF5[lenF5]; }
  void addvflen(const int32_t flen) { ++vflen[flen]; }

  void outputDistFile(const std::string &prefix, const uint64_t nread) {
    std::string outputfile = prefix + ".ReadLengthDist.csv";
    std::ofstream out(outputfile);
    printVector(out, vlenF3, "F3", nread);
    if(pairedend) {
      out << std::endl;
      printVector(out, vlenF5, "F5", nread);
    }
    out.close();
    
    if(pairedend) {
      outputfile = prefix + ".FragmentLengthDist.csv";
      std::ofstream out(outputfile);
      printVector(out, vflen, "Fragmemt", nread);
    }
  }
};

class SeqWigStats: public SeqStats {
 public:
  WigStats ws;
  
 SeqWigStats(std::string &s, int32_t l=0, int32_t binsize=0): SeqStats(s, l, binsize) {}
};

class SeqStatsGenome {
  std::string name;
  double depth;
  double sizefactor;

  std::vector<bed> vbed;
  void readGenomeTable(const std::string &gt, const int binsize);

 public:
  std::vector<SeqWigStats> chr;
  std::vector<MyMthread::chrrange> vsepchr;
  FragmentLengthDist dflen;
  
 SeqStatsGenome(const MyOpt::Variables &values):
  name("Genome"), depth(0), sizefactor(0), dflen(values)
  {
    readGenomeTable(values["gt"].as<std::string>(), values["binsize"].as<int32_t>());
    if(values.count("mptable")) {
      for(auto &x: chr) x.getMptable(values["mptable"].as<std::string>());
    }
    
    // Greekchr
    for(auto &x: chr) {
      if(x.getname() == "I") {
	for(auto &x:chr) x.Greekchron();
	break;
      }
    }

    // sepchr
    vsepchr = MyMthread::getVsepchr(getlen(), chr, values["threads"].as<int32_t>());
    
#ifdef DEBUG
    std::cout << "chr\tautosome" << std::endl;
    for(auto &x: chr) std::cout << x.getname() << "\t" << x.isautosome() << std::endl;
    for(uint32_t i=0; i<vsepchr.size(); i++)
      std::cout << "thread " << (i+1) << ": " << vsepchr[i].s << "-" << vsepchr[i].e << std::endl;
    printReadstats();
#endif
  }

  std::string getname() const { return name; }
  uint64_t getlen() const {
    uint64_t len(0);
    for(auto &x:chr) len += x.getlen();
    return len;
  }
  uint64_t getlenmpbl() const {
    uint64_t len_mpbl(0);
    for(auto &x:chr) len_mpbl += x.getlenmpbl();
    return len_mpbl;
  }
  double getpmpbl() const {
    return getratio(getlenmpbl(), getlen());
  }
  int32_t getnbin() const {
    int32_t nbin(0);
    for(auto &x: chr) nbin += x.getnbin();
    return nbin;
  }
  uint64_t getnbp() const {
    uint64_t nbp(0);
    for(auto &x: chr) nbp += x.getnbp();
    return nbp;
  }
  uint64_t getncov() const {
    uint64_t ncov(0);
    for(auto &x: chr) ncov += x.getncov();
    return ncov;
  }
  uint64_t getncovnorm() const {
    uint64_t ncovnorm(0);
    for(auto &x: chr) ncovnorm += x.getncovnorm();
    return ncovnorm;
  }
  
  uint64_t getnread (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread(strand);
    return nread;
  }
  uint64_t getnread_nonred (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_nonred(strand);
    return nread;
  }
  uint64_t getnread_red (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_red(strand);
    return nread;
  }
  uint64_t getnread_rpm (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_rpm(strand);
    return nread;
  }
  uint64_t getnread_afterGC (const Strand::Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_afterGC(strand);
    return nread;
  }
  uint64_t getnread_inbed() const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_inbed();
    return nread;
  }
  double getFRiP() const {
    return getratio(getnread_inbed(), getnread_nonred(Strand::BOTH));
  }
  void setdepth(const double d) { depth = d; }
  double getdepth() const { return depth; }
  double getsizefactor()const { return sizefactor; }

  void setsizefactor(const double w) { sizefactor = w; }
  void setbed(const std::string &bedfilename) {
    isFile(bedfilename);
    vbed = parseBed<bed>(bedfilename);
    //    printBed(vbed);
  }
  const std::vector<bed> & getvbedref() const { return vbed; }
  
  void setF5All(const int32_t flen) {
    for (auto &x:chr) x.setF5(flen);
  }

  void printReadstats() const {
    std::cout << "name\tlength\tlen_mpbl\tread num\tnonred num\tred num\tnormed\tafterGC\tdepth" << std::endl;
    printSeqStats(*this);
    for(auto &x: chr) printSeqStats(x);
  }
};

#endif /* _MAPFILE_HPP_ */
