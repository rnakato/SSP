/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _BEDFORMAT_HPP_
#define _BEDFORMAT_HPP_

#include "util.hpp"

namespace NAKATO {
class bed {
 public:
  std::string chr;
  int32_t start;
  int32_t end;
  int32_t summit;
 bed(): start(0), end(0), summit(0) {}
  virtual ~bed(){}
  bed(const std::string &c, const int32_t s, const int32_t e, const int32_t _summit=0):
    chr(rmchr(c)), start(s), end(e)
  {
    if (_summit) summit = _summit;
    else summit = (start + end)/2;
  }
  explicit bed(const std::vector<std::string> &s) {
    if(s.size() < 3) {
      std::cerr << "\nWarning: Bed site size < 3." << std::endl;
      return;
    }
    chr = rmchr(s[0]);
    start = stoi(s[1]);
    end = stoi(s[2]);
    summit = (start + end)/2;
  }
  void print() const { std::cout << "chr" << chr << "\t" << start  << "\t" << end ; }
  void printHead() const { std::cout << "chromosome\tstart\tend"; }
  int32_t length() const { return abs(end - start); }
  std::string getSiteStr() const {
    return "chr" + chr + "-" + std::to_string(start) + "-" + std::to_string(end);
  }

};
}

#endif  // _BEDFORMAT_HPP_
