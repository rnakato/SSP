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

 public:
  chrsize(const std::string &n, const int32_t l):
    refname(n), name(rmchr(n)), len(l)
  {}

  const std::string & getrefname() const { return refname; }
  const std::string & getname() const { return name; }
  int32_t getlen() const { return len; }
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
 var(): name(""), val(0), limlow(0), limup(0), isupper(false){}
  var(std::string &str, T low):       name(str), val(0), limlow(low), limup(0), isupper(false) {}
  var(std::string &str, T low, T up): name(str), val(0), limlow(low), limup(up), isupper(true) {}
  void set(T n) {
    if (isupper && (n<limlow || n>limup)) {
      std::cout << "Error : variable " << name << " should be " << limlow << "<= and <=" << limup << "." << std::endl;
    }else if (!isupper && n<limlow ) {
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
    std::cout << name << "\t"
	      << len << "\t"
	      << nbin << "\t"
	      << len_mpbl << "\t"
	      << p_mpbl << "\t"
	      << gcov << std::endl;
  }
};

class Fragment {
public:
  std::string chr;
  int32_t F3;
  Strand::Strand strand;
  int32_t fraglen;
  int32_t readlen_F3;

  Fragment(): chr(""), F3(0), strand(Strand::FWD), fraglen(0), readlen_F3(0) {}

  void addSAM(const std::string &_chrname, const int32_t _readlen, const int32_t position,
	      const int32_t isize, const bool _strand, const bool ispaired)
  {
    chr = rmchr(_chrname);
    readlen_F3 = _readlen;
    if (ispaired) fraglen = abs(isize);
    if (_strand) {  // 0: forward 1: reverse
      strand = Strand::REV;
      F3 = position + readlen_F3;
    } else {
      strand = Strand::FWD;
      F3 = position;
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
    if (frag.strand == Strand::FWD) F5 = frag.F3 + frag.fraglen;
    else                            F5 = frag.F3 - frag.fraglen;
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
