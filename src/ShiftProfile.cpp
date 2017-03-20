/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <map>
#include <numeric>
#include <boost/thread.hpp>
#include "ShiftProfile.hpp"
#include "ShiftProfile_p.hpp"
#include "../common/inline.hpp"
#include "../common/statistics.hpp"

void addmp(std::map<int32_t, double> &mpto, const std::map<int32_t, double> &mpfrom, double w)
{
  for(auto x: mpfrom) {
    if(!std::isnan(x.second)) mpto[x.first] += x.second * w;
  }
}

double getJaccard(int32_t step, int32_t to, int32_t xysum, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  int32_t xy(0);
  for(int32_t j=mp_from; j<to; ++j) if(fwd[j] * rev[j+step]) xy += std::max(fwd[j], rev[j+step]);
  return getratio(xy, xysum-xy);
}

void genThreadJacVec(ReadShiftProfile &chr, int32_t ng_to, int32_t xysum, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev, int32_t s, int32_t e, boost::mutex &mtx)
{
  for(int32_t step=s; step<e; ++step) {
    chr.setmp(step, getJaccard(step, chr.width-ng_to, xysum, fwd, rev), mtx);
  }
}

void shiftJacVec::setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  int32_t xx = accumulate(fwd.begin(), fwd.end(), 0);
  int32_t yy = accumulate(rev.begin(), rev.end(), 0);

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint32_t i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(&genThreadJacVec, boost::ref(chr), ng_to, xx+yy, boost::cref(fwd), boost::cref(rev), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for(int32_t step=ng_from; step<ng_to; step+=ng_step) {
    chr.nc[step] = getJaccard(step, chr.width-ng_to, xx+yy, fwd, rev);
  }
}

void genThreadCcp(ReadShiftProfile &chr, int32_t ng_to, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev, double mx, double my, const int32_t s, const int32_t e, boost::mutex &mtx)
{
  for(int32_t step=s; step<e; ++step) {
    double xy(0);
    for(int32_t j=mp_from; j<chr.width-ng_to; ++j) {
      xy += (fwd[j] - mx) * (rev[j+step] - my);
    }
    chr.setmp(step, xy, mtx);
  }
}

void shiftCcp::setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  if(chr.width < ng_to) {
    std::cerr << "\nerror: chromosome length " << chr.width << " is shorter than background distance " << ng_to << std::endl;
    std::cerr << "please specify shorter length with --ng_from and --ng_to options." << std::endl;
    exit(0);
  }
  
  MyStatistics::moment<int8_t> x(fwd, mp_from, chr.width - ng_to);
  MyStatistics::moment<int8_t> y(rev, mp_from, chr.width - ng_to);

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint32_t i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(genThreadCcp, boost::ref(chr), ng_to, boost::cref(fwd), boost::cref(rev), x.getmean(), y.getmean(), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for(int32_t step=ng_from; step<ng_to; step+=ng_step) {
    double xy(0);
    for(int32_t j=mp_from; j<chr.width-ng_to; ++j) xy += (fwd[j] - x.getmean()) * (rev[j+step] - y.getmean());
    chr.nc[step] = xy;
  }

  double val = 1/(x.getsd() * y.getsd() * (chr.width - ng_to - mp_from - 1));
  std::cout<< val << std::endl;
  for(auto &x: mp) x.second *= val;
  for(auto &x: nc) x.second *= val;
}

void shiftJacBit::setDist(ReadShiftProfile &chr, boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  /*

  rev >>= 150; //getlenF3();
  std::vector<int8_t> array(fwd.size()/100 +1, 0);
  for(uint32_t i=0; i<fwd.size()-getlenF3(); ++i) {
    if(fwd.test(i) && rev.test(i)) {
      //      fwd.reset(i);
      //rev.reset(i);
      ++array[i/100]; //std::cout << i << std::endl;
    }
  }
  
  std::cerr << (int)MyStatistics::getPercentile(array, 0.95) << "\t"<< (int)MyStatistics::getPercentile(array, 0.98) << std::endl;
  int max = MyStatistics::getPercentile(array, 0.98);
  for(uint32_t i=0; i<fwd.size()/100; ++i) {
    if(array[i]>=2) {
      std::cout << (i*100) << "\t"<< ((int)array[i]) << std::endl;
      for(int32_t j=i*100; j<(i+1)*100; ++j) {
	fwd.reset(j);
	rev.reset(j);
      }
    }
  }
  rev <<= getlenF3();
  exit(0);*/
  
  double xysum(fwd.count() + rev.count());
  
  rev <<= mp_from;
  for(int32_t step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    int32_t xy((fwd & rev).count());
    chr.mp[step] = xy/(xysum-xy);
  }
  
  rev >>= (ng_from - mp_to);

  for(int32_t step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    int32_t xy((fwd & rev).count());
    chr.nc[step] = xy/(xysum-xy);
  }
}

void shiftHamming::setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  rev <<= mp_from;
  for(int32_t step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    chr.mp[step] = (fwd ^ rev).count();
  }
  rev >>= (ng_from - mp_to);
  
  for(int32_t step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    chr.nc[step] = (fwd ^ rev).count();
  }
}

std::vector<int8_t> genVector(const std::vector<Read> &vReadref, const int32_t start, const int32_t end)
{
  std::vector<int8_t> array(end-start, 0);
  for (auto &x: vReadref) {
    if(!x.duplicate && my_range(x.F3, start, end-1)) ++array[x.F3 - start];
  }
  return array;
}

boost::dynamic_bitset<> genBitset(const std::vector<Read> &vReadref, const int32_t start, const int32_t end)
{
  boost::dynamic_bitset<> array(end-start);
  for (auto &x: vReadref) {
    if(!x.duplicate && my_range(x.F3, start, end-1))
      array.set(x.F3 - start);
  }
  return array;
}

namespace {
  void setSSPstats(SSPstats &p, const double bu, const double nsc, const double rsc)
  {
    p.setnsc(nsc);
    p.setrsc(rsc);
    p.setbu(bu);
  }
}

template <class T>
void genThread(T &dist, const SeqStatsGenome &genome, uint32_t chr_s, uint32_t chr_e, const std::string &prefix, const bool output_eachchr)
{
  for(uint32_t i=chr_s; i<=chr_e; ++i) {
    std::cout << genome.chr[i].getname() << ".." << std::flush;

    dist.execchr(genome, i);
    dist.chr[i].setflen(dist.name);
    
    std::string filename = prefix + "." + genome.chr[i].getname() + ".csv";
    if(output_eachchr) dist.outputmpChr(filename, i);
    //    exit(0);
  }
}

template <class T>
void makeProfile(SSPstats &sspst, SeqStatsGenome &genome, const std::string &head, const std::string &typestr)
{
  DEBUGprint("makeProfile: " + typestr );
  T dist(sspst, genome);
  dist.printStartMessage();
  
  boost::thread_group agroup;
  boost::mutex mtx;

  std::string prefix(head + "." + typestr);
  if(typestr == "hdp" || typestr == "jaccard") {
    for(size_t i=0; i<genome.vsepchr.size(); i++) {
      agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(genome), genome.vsepchr[i].s, genome.vsepchr[i].e, boost::cref(prefix), sspst.isEachchr()));
    }
    agroup.join_all();
  } else {
    genThread(dist, genome, 0, genome.chr.size()-1, prefix, sspst.isEachchr());
  }

  for(size_t i=0; i<genome.chr.size(); ++i) {
    if(genome.chr[i].isautosome()) dist.addmp2genome(i);
  }

  dist.setflen(dist.name);
  genome.dflen.setflen_ssp(dist.getnsci());

  std::string prefix2 = head + "." + typestr;
  dist.outputmpGenome(prefix2);

  if(typestr == "jaccard") setSSPstats(sspst, dist.getbackgroundUniformity(), dist.getnsc(), dist.getrsc());

  DEBUGprint("makeProfile: " + typestr + " done.");
  return;
}

void strShiftProfile(SSPstats &sspst, SeqStatsGenome &genome, const std::string &head, const std::string &typestr)
{
  DEBUGprint("strShiftProfile...");
  
  if(typestr=="exjaccard")    makeProfile<shiftJacVec>(sspst, genome, head, typestr);
  else if(typestr=="jaccard") makeProfile<shiftJacBit>(sspst, genome, head, typestr);
  else if(typestr=="ccp")     makeProfile<shiftCcp>(sspst, genome, head, typestr);
  else if(typestr=="hdp")     makeProfile<shiftHamming>(sspst, genome, head, typestr);

  DEBUGprint("strShiftProfile done.");
  return;
}
