/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _WIGSTATS_H_
#define _WIGSTATS_H_

#include <vector>
#include <fstream>
#include "macro.h"

class WigArray {
  std::vector<int32_t> array;
  double geta;
  
  WigArray(){}

  template <class T>
    double rmGeta(T val) const { return val/geta; } 
  template <class T>
    double addGeta(T val) const { return val*geta; }
  void checki(const size_t i) const {
    if(i>=size()) PRINTERR("Invalid i for WigArray: " << i << " > " << size());
  }

 public:
 WigArray(const size_t num, const int32_t val): array(num, val), geta(1000.0) {}
  ~WigArray(){}

  size_t size() const { return array.size(); }
  double getval(const size_t i) const {
    checki(i);
    return rmGeta(array[i]);
  }
  void setval(const size_t i, const double val) {
    checki(i);
    array[i] = addGeta(val);
  }
  void addval(const size_t i, const double val) {
    checki(i);
    array[i] += addGeta(val);
  }
  void multipleval(const size_t i, const double val) {
    checki(i);
    array[i] *= val;
  }

  double getPercentile(double per) const {
    int32_t v95(MyStatistics::getPercentile(array, per));
    return rmGeta(v95);
  }

  void outputAsWig(std::ofstream &out, const int32_t binsize) const {
    for(size_t i=0; i<array.size(); ++i) {
      if(array[i]) out << boost::format("%1%\t%2%\n") % (i*binsize+1) % rmGeta(array[i]);
    }
  }
  void outputAsBedGraph(std::ofstream &out, const int32_t binsize, const std::string &name, const uint64_t chrend) const {
    uint64_t e;
    for(size_t i=0; i<array.size(); ++i) {
      if(i==array.size() -1) e = chrend; else e = (i+1) * binsize;
      if(array[i]) out << boost::format("%1% %2% %3% %4%\n") % name % (i*binsize) % e % rmGeta(array[i]);
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

#endif /* _MAPFILECLASS_H_ */
