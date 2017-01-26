/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <fstream>
#include <numeric>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include "mapfileclass.h"
#include "util.h"
#include "readdata.h"
#include "mthread.h"

namespace MyOpt {
  using Variables = boost::program_options::variables_map;
  using Opts      = boost::program_options::options_description;
}

void printDist(std::ofstream &out, const std::vector<int32_t> v, const std::string str, const uint64_t nread);

class SeqStatsGenome {
  std::string name;
  double depth;
  double sizefactor;

  std::vector<bed> vbed;
  void readGenomeTable(const std::string &gt, const int binsize) {
    std::vector<std::string> v;
    std::string lineStr;
    std::ifstream in(gt);
    if(!in) PRINTERR("Could nome open " << gt << ".");
    
    while (!in.eof()) {
      getline(in, lineStr);
      if(lineStr.empty() || lineStr[0] == '#') continue;
      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      SeqStats s(v[0], stoi(v[1]), binsize);
      chr.push_back(s);
    }
    return;
  }

 public:
  WigStats ws;
  std::vector<SeqStats> chr;
  std::vector<MyMthread::chrrange> vsepchr;
  
 SeqStatsGenome(const MyOpt::Variables &values): name("Genome"), depth(0), sizefactor(0) {
    readGenomeTable(values["gt"].as<std::string>(), values["binsize"].as<int32_t>());
    if(values.count("mp")) getMpbl(values["mp"].as<std::string>(), chr);
    else if(values.count("mptable")) getMpbltable(values["mptable"].as<std::string>(), chr);

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
  
  uint64_t getnread (const Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread(strand);
    return nread;
  }
  uint64_t getnread_nonred (const Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_nonred(strand);
    return nread;
  }
  uint64_t getnread_red (const Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_red(strand);
    return nread;
  }
  uint64_t getnread_rpm (const Strand strand) const {
    uint64_t nread(0);
    for(auto &x:chr) nread += x.getnread_rpm(strand);
    return nread;
  }
  uint64_t getnread_afterGC (const Strand strand) const {
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
    return getratio(getnread_inbed(), getnread_nonred(STRAND_BOTH));
  }
  void setdepth(const double d) { depth = d; }
  double getdepth() const { return depth; }
  double getsizefactor()const { return sizefactor; }

  void setsizefactor(const double w) { sizefactor = w; }
  void setbed(const std::string bedfilename) {
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
  Greekchr(false),
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
    out << samplename << "\t" << genome.getnread(STRAND_BOTH) << "\t" << genome.getnread_nonred(STRAND_BOTH) << "\t"
	<< lenF3 << "\t" << eflen << "\t";
    sspst.print(out);
  }
  
  void setmaxGC(const int32_t m) { maxGC = m; }
  int32_t getmaxGC() const {return maxGC; }

  void lackOfRead4Complexity_on() { lackOfRead4Complexity = true; }
  void lackOfRead4GenomeCov_on() { lackOfRead4GenomeCov = true; }
  bool islackOfRead4GenomeCov() const { return lackOfRead4GenomeCov; };
  void setthre4filtering(const MyOpt::Variables &values) {
    if(values["thre_pb"].as<int32_t>()) thre4filtering = values["thre_pb"].as<int32_t>();
    else {
      int32_t thre = getratio(genome.getnread(STRAND_BOTH), genome.getlenmpbl()) * 10;
      thre4filtering = std::max(1, thre);
    }
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
      if(x.getname() == frag.chr) {
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
    printDist(out, vlenF3, "F3", genome.getnread(STRAND_BOTH));
    if(values.count("pair")) printDist(out, vlenF5, "F5", genome.getnread(STRAND_BOTH));
    out.close();
    
    if(values.count("pair")) {
      outputfile = oprefix + ".fraglen_dist.xls";
      std::ofstream out(outputfile);
      printDist(out, vflen, "Fragmemt", genome.getnread(STRAND_BOTH));
    }
  }

  void printComplexity(std::ofstream &out) const {
    if(lackOfRead4Complexity) out << boost::format("Library complexity: (%1$.3f) (%2%/%3%)\n") % complexity() % nt_nonred % nt_all;
    else out << boost::format("Library complexity: %1$.3f (%2%/%3%)\n") % complexity() % nt_nonred % nt_all;
  }
  double complexity() const { return getratio(nt_nonred, nt_all); }

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
