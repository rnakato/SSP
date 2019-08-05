/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SEQ_H_
#define _SEQ_H_

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

std::string rmchr(const std::string &chr);

namespace Strand {
  enum Strand {FWD, REV, BOTH};
}

class chrsize {
  std::string refname;
  std::string name;
  int32_t len;
  //  bool Greekchr;

 public:
  chrsize(const std::string &n, const int32_t l):
    refname(n), name(rmchr(n)), len(l) //, Greekchr(false)
  {}

  const std::string & getrefname() const { return refname; }
  const std::string & getname() const { return name; }
  /*  const std::string & getrefname_GreekToInt() const {
    std::string chrnum;
    if(name=="chrI")         chrnum = "chr1";
    else if(name=="chrII")   chrnum = "chr2";
    else if(name=="chrIII")  chrnum = "chr3";
    else if(name=="chrIV")   chrnum = "chr4";
    else if(name=="chrV")    chrnum = "chr5";
    else if(name=="chrVI")   chrnum = "chr6";
    else if(name=="chrVII")  chrnum = "chr7";
    else if(name=="chrVIII") chrnum = "chr8";
    else if(name=="chrIX")   chrnum = "chr9";
    else if(name=="chrX")    chrnum = "chr10";
    else if(name=="chrXI")   chrnum = "chr11";
    else if(name=="chrXII")  chrnum = "chr12";
    else if(name=="chrXIII") chrnum = "chr13";
    else if(name=="chrXIV")  chrnum = "chr14";
    else if(name=="chrXV")   chrnum = "chr15";
    else if(name=="XVIchr")  chrnum = "chr16";
    return chrnum;
    }*/
  /*  const std::string getname_GreekToInt(const std::string &name) const {
    std::string chrnum;
    if(name=="I")         chrnum = "1";
    else if(name=="II")   chrnum = "2";
    else if(name=="III")  chrnum = "3";
    else if(name=="IV")   chrnum = "4";
    else if(name=="V")    chrnum = "5";
    else if(name=="VI")   chrnum = "6";
    else if(name=="VII")  chrnum = "7";
    else if(name=="VIII") chrnum = "8";
    else if(name=="IX")   chrnum = "9";
    else if(name=="X")    chrnum = "10";
    else if(name=="XI")   chrnum = "11";
    else if(name=="XII")  chrnum = "12";
    else if(name=="XIII") chrnum = "13";
    else if(name=="XIV")  chrnum = "14";
    else if(name=="XV")   chrnum = "15";
    else if(name=="XVI")  chrnum = "16";
    return chrnum;
    }*/
  int32_t getlen() const { return len; }
  /*  void Greekchron() {
    name = getname_GreekToInt(name);
    //Greekchr = true;
    }*/
};

class range {
 public:
  int32_t start;
  int32_t end;
 range(): start(0), end(0) {}
 range(int32_t s, int32_t e): start(s), end(e) {}
  int32_t getlen() const { return end-start; }
};

template <class T>
class var {
  std::string name;
  T val;
  T limlow;
  T limup;
  bool isupper;
 public:
 var(): name(""), val(0), limlow(0), limup(0){}
  var(std::string &str, T low):       name(str), val(0), limlow(low), limup(0), isupper(false) {}
  var(std::string &str, T low, T up): name(str), val(0), limlow(low), limup(up), isupper(true) {}
  void set(T n) {
    if(isupper && (n<limlow || n>limup)) {
      std::cout << "Error : variable " << name << " should be " << limlow << "<= and <=" << limup << "." << std::endl;
    }else if(!isupper && n<limlow ) {
      std::cout << "Error : variable " << name << " should be >=" << limlow << "." << std::endl;
    }
    else val=n;
  }
  operator T() const { return val; }
};

class fasta {
 public:
  std::string name;
  uint64_t len, len_mpbl;
  int32_t nbin;
  double p_mpbl;  /* mappability */
  double gcov;    /* genome coverage for bin */
  fasta (std::string &str, int32_t l=0): name(str), len(l), len_mpbl(0), nbin(0), p_mpbl(0), gcov(0) {}
  explicit fasta (std::vector<std::string> &v): name(v[0]), len(stoi(v[1])), len_mpbl(0), nbin(0), p_mpbl(0), gcov(0) {}
  void print() const {
    std::cout << name << "\t" << len << "\t" << nbin << "\t" << len_mpbl << "\t"<< p_mpbl << "\t" << gcov << std::endl;
  }
};

class Fragment {
public:
  std::string chr;
  int32_t F3;
  Strand::Strand strand;
  int32_t fraglen;
  int32_t readlen_F3;

  Fragment(): F3(0), fraglen(0), readlen_F3(0) {}
  void addSAM(const std::vector<std::string> &v, const bool pair, const int32_t sv) {
   chr = rmchr(v[2]);
   readlen_F3 = v[9].length();
   if(pair) fraglen = abs(stoi(v[8]));
   if(sv&16) {
     strand = Strand::REV;
     F3 = stoi(v[3]) + readlen_F3 -1;
   } else {
     strand = Strand::FWD;
     F3 = stoi(v[3]) -1;
   }
 }
 void print() const {
#ifdef PRINTFRAGMENT
   std::cout << "chr:"       << chr
	     << "\tposi:"    << F3
	     << "\tstrand:"  << strand
	     << "\tfraglen:" << fraglen
	     <<"\treadlen:"  << readlen_F3
	     << std::endl;
#endif
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

  Read(const Fragment &frag, const int32_t len):
    weight(WeightNum), F3(frag.F3), duplicate(0), inpeak(0)
  {
    if(frag.strand == Strand::FWD) F5 = frag.F3 + frag.fraglen;
    else                           F5 = frag.F3 - frag.fraglen;
    F3 = std::max(F3, 0);
    F5 = std::max(F5, 0);
    F3 = std::min(F3, len-1);
    F5 = std::min(F5, len-1);
  }
  double getWeight() const {
    return weight/static_cast<double>(WeightNum);
  }
  void multiplyWeight(const double w) { weight *= w; }
  void print() const {
#ifdef PRINTREAD
    std::cout << "F3:"      << F3
	      << "\tF5:"    << F5
	      << "\tweight:"<< weight
	      << "\tduplicate:" << duplicate
	      << "\tinpeak:"    << inpeak
	      << std::endl;
#endif
  }
};


#endif  // _SEQ_H_
