/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _BEDFORMAT_HPP_
#define _BEDFORMAT_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <boost/algorithm/string.hpp>

std::string rmchr(const std::string &chr);
bool isStr(std::string, std::string);

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
  bed(const std::vector<std::string> &s):
    chr(rmchr(s[0])), start(stoi(s[1])), end(stoi(s[2])), summit((start + end)/2)
  {}
  void print() const { std::cout << "chr" << chr << "\t" << start  << "\t" << end ; }
  void printHead() const { std::cout << "chromosome\tstart\tend"; }
  int32_t length() const { return abs(end - start); }
};

class bed12 : public bed {
 public:
  std::string name;
  int32_t score;
  std::string strand;
  int32_t thickStart;
  int32_t thickEnd;
  std::string itemRgb;
  int32_t blockCount;
  int32_t blockSizes;
  int32_t blockStarts;
 bed12(): bed() {}
 bed12(std::vector<std::string> s): bed(s) {
   int32_t num = s.size();
   if(num > 3)  name        = s[3];
   if(num > 4)  score       = stoi(s[4]);
   if(num > 5)  strand      = s[5];
   if(num > 6)  thickStart  = stoi(s[6]);
   if(num > 7)  thickEnd    = stoi(s[7]);
   if(num > 8)  itemRgb     = s[8];
   if(num > 9)  blockCount  = stoi(s[9]);
   if(num > 10) blockSizes  = stoi(s[10]);
   if(num > 11) blockStarts = stoi(s[11]);
 }
 void print() const {
   std::cout << chr << "\t" << start << "\t" << end << "\t"
	<< name << "\t" << score << "\t" << strand << "\t"
	<< thickStart << "\t" << thickEnd << "\t" << itemRgb << "\t"
	<< blockCount << "\t" << blockSizes << "\t" << blockStarts;
 }
 void printHead() const {
   std::cout << "chromosome\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
 }
};

class macsxls : public bed {
 public:
  int32_t len;
  int32_t summit;
  double pileup;
  double p;
  double enrich;
  double q;
  std::string name;

 macsxls(): bed() {}
 macsxls(std::vector<std::string> s): bed(s) {
   len    = stoi(s[3]);
   summit = stoi(s[4]);
   pileup = stod(s[5]);
   p      = stod(s[6]);
   enrich = stod(s[7]);
   q      = stod(s[8]);
   name   = s[9];
 }
 void print() const {
   std::cout << chr << "\t" << start << "\t" << end << "\t"
	<< len << "\t" << summit << "\t" << pileup << "\t"
	<< p << "\t" << enrich << "\t" << q << "\t" << name;
 }
 void printHead () const {
   std::cout << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname";
 }
};

class Peak : public bed {
 public:
  int32_t summit;
  double pileup;
  double enrich;
  double p_inter, p_enr;
  double q;

  Peak(){}
  Peak(const std::string &c, const int32_t s, const int32_t e, const double val, const double p):
    bed(c,s,e), summit(s), pileup(val), enrich(0), p_inter(p), p_enr(0), q(0) {}
  void renew(const int32_t i, const double val, const double p) {
    end = i;
    pileup += val;
    if(p_inter > p) {
      p_inter = p;
      summit = i;
    }
  }
  void print(std::ofstream &out, const int32_t id, const int32_t binsize) const {
    out << chr << "\t" << start*binsize << "\t" << end*binsize << "\t"
	<< ((end - start +1)*binsize-1) << "\t"
	<< (summit*binsize -binsize/2) << "\t" << pileup << "\t"
	<< p_inter << "\t" << enrich << "\t" << q << "\tpeak " << id << std::endl;
  }
  void printHead (std::ofstream &out) const {
    out << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname" << std::endl;
  }
};

template <class T>
std::vector<T> parseBed(const std::string &fileName)
{
  std::vector<T> vbed;
  std::ifstream in(fileName);
  if(!in) {
    std::cerr << "Error: BED file " << fileName << " does not exist." << std::endl;
    std::exit(1);
  }

  while (!in.eof()) {
    std::string lineStr;
    getline(in, lineStr);

    if(lineStr.empty() || lineStr[0] == '#') continue;
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v[1] == "start") continue;
    vbed.emplace_back(v);
  }

  return vbed;
}

template <class T>
void printBed(const std::vector<T> &vbed)
{
  for (auto &x: vbed) {
    x.print();
    std::cout << std::endl;
  }
  std::cout << "Total number: " << vbed.size() << std::endl;
  return;
}

class Interaction {
  double val;
 public:
  bed first;
  bed second;
  Interaction(): val(0) {}
  Interaction(const std::string &c1, const int32_t s1, const int32_t e1,
	      const std::string &c2, const int32_t s2, const int32_t e2,
	      const double v=0):
    val(v), first(c1, s1, e1), second(c2, s2, e2)
  {}
  Interaction(const bed &b1, const bed &b2, const double _val):
    val(_val), first(b1), second(b2)
  {}
  
  double getval() const { return val; }
  void print() const {
    first.print();
    std::cout << "\t";
    second.print();
    std::cout << "\t" << val << std::endl;
  }
};

class InteractionSet {
  std::vector<Interaction> vinter;
  double maxval;
  std::string label;

  void setAsMango(const std::string &lineStr) {
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v.size() < 8) {
      std::cerr << "Warning: " << lineStr << " does not contain 8 columns." << std::endl;
      return;
    }
    try {
      double p(1e-12);
      if(stod(v[7])) p = stod(v[7]);
      double val(-log10(p));
      vinter.emplace_back(bed({v[0], v[1], v[2]}),
			  bed({v[3], v[4], v[5]}),
			  val);
      maxval = std::max(val, maxval);
    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      exit(0);
    }
  }
  void setAsHICCUPS(const std::string &lineStr) {
    if(isStr(lineStr, "color")) return;
    std::vector<std::string> v;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v.size() < 19) {
      std::cerr << "Warning: " << lineStr << " does not contain 8 columns." << std::endl;
      return;
    }

    try {
      double val(-log10(stod(v[13])));
      //      std::cout << val << "\t" << stod(v[13]) << "\t" << v[13]  << std::endl;
      vinter.emplace_back(bed(v[0], stoi(v[1]), stoi(v[2]), stoi(v[17])),
			  bed(v[3], stoi(v[4]), stoi(v[5]), stoi(v[18])),
			  val);
      maxval = std::max(val, maxval);
    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      exit(0);
    }
  }
  
public:
  InteractionSet(const std::string &fileName, const std::string &l, const std::string &tool):
    maxval(0), label(l)
  {
    std::ifstream in(fileName);
    if(!in) {
      std::cerr << "Error: Interaction file " << fileName << " does not exist." << std::endl;
      std::exit(1);
    }

    while (!in.eof()) {
      std::string lineStr;
      getline(in, lineStr);
      if(lineStr.empty() || lineStr[0] == '#') continue;
      if (tool == "mango") setAsMango(lineStr);
      else setAsHICCUPS(lineStr);
    }

    print();
  }
  const std::vector<Interaction> & getvinter() const { return vinter; }
  const std::string & getlabel() const { return label; }
  double getmaxval() const { return maxval; }
  size_t getnum() const { return vinter.size(); }
  void print() const {
    printBed(vinter);
    std::cout << "maxval: " << maxval << std::endl;
  }
};

template <class T>
std::unordered_map<std::string, std::vector<T>> parseBed_Hash(const std::string &fileName)
{
  std::unordered_map<std::string, std::vector<T>> bedmap;
  std::ifstream in(fileName);
  if(!in) {
    std::cerr << "Error: BED file does not exist." << std::endl;
    std::exit(1);
  }

  std::string lineStr;
  std::vector<std::string> v;
  while (!in.eof()) {
    getline(in, lineStr);

    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v[1] == "start") continue;
    bedmap[rmchr(v[0])].emplace_back(bed(v));
  }

  return bedmap;
}

template <class T>
void printBed_Hash(const std::unordered_map<std::string, std::vector<T>> &mp)
{
  for(auto vbed: mp) printBed(vbed.second);
  return;
}


#endif  // _BEDFORMAT_HPP_
