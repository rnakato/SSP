/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <map>
#include <boost/thread.hpp>
#include "ShiftProfile.hpp"
#include "ShiftProfile_p.hpp"
#include "../common/inline.hpp"
#include "../common/statistics.hpp"

void addmp(std::map<int32_t, double> &mpto, const std::map<int32_t, double> &mpfrom, double w)
{
  for (auto x: mpfrom) {
    if (!std::isnan(x.second)) mpto[x.first] += x.second * w;
  }
}

double shiftJacVec::getJaccard(const int32_t step, const int32_t to,
			       const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  int64_t xy(0);
  for (int32_t j=mp_from; j<to; ++j) {
    if (std::min(fwd[j], rev[j+step])) xy += std::max(fwd[j], rev[j+step]);
  }
  return xy;
}

void shiftJacVec::genThreadJacVec(ReadShiftProfile &chr,
				  const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev,
				  const int32_t s, const int32_t e, boost::mutex &mtx)
{
  for (int32_t step=s; step<e; ++step) {
    chr.setmp(step, getJaccard(step, chr.width-ng_to, fwd, rev), mtx);
  }
}

void shiftJacVec::setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  int64_t bothsum(0);
  for (auto x: fwd) bothsum += x;
  for (auto x: rev) bothsum += x;

  boost::thread_group agroup;
  boost::mutex mtx;
  for (uint32_t i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(&shiftJacVec::genThreadJacVec, this, boost::ref(chr), boost::cref(fwd), boost::cref(rev), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for (int32_t step=ng_from; step<ng_to; step+=ng_step) {
    chr.nc[step] = getJaccard(step, chr.width-ng_to, fwd, rev);
  }
}

void shiftCcp::genThread(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev,
			 const double mx, const double my, const int32_t s, const int32_t e, boost::mutex &mtx)
{
  for (int32_t step=s; step<e; ++step) {
    double xy(0);
    for (int32_t j=mp_from; j<chr.width-ng_to; ++j) {
      xy += (fwd[j] - mx) * (rev[j+step] - my);
    }
    chr.setmp(step, xy, mtx);
  }
}

void shiftCcp::setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  MyStatistics::moment<int8_t> x(fwd, mp_from, chr.width - ng_to);
  MyStatistics::moment<int8_t> y(rev, mp_from, chr.width - ng_to);

  boost::thread_group agroup;
  boost::mutex mtx;
  for (uint32_t i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(&shiftCcp::genThread, this, boost::ref(chr), boost::cref(fwd), boost::cref(rev), x.getmean(), y.getmean(), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for (int32_t step=ng_from; step<ng_to; step+=ng_step) {
    double xy(0);
    for (int32_t j=mp_from; j<chr.width-ng_to; ++j) xy += (fwd[j] - x.getmean()) * (rev[j+step] - y.getmean());
    //    printf("step %d, xy %f\n", step, xy);
    chr.nc[step] = xy;
  }

  double val = 1/(x.getsd() * y.getsd() * (chr.width - ng_to - mp_from - 1));
  for (auto &x: chr.mp) x.second *= val;
  for (auto &x: chr.nc) x.second = std::abs(x.second * val); //*= val;
}

void shiftJacBit::setDist(ReadShiftProfile &chr, boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  double xysum(fwd.count() + rev.count());

  clock_t t1,t2;
  t1 = clock();

  rev <<= mp_from;
  for (int32_t step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    int32_t xy((fwd & rev).count());
    chr.mp[step] = xy/(xysum-xy);
    //    printf("step %d %f %d %f\n", step, chr.mp[step], xy, xysum);
  }
  t2 = clock();
  PrintTime(t1, t2, "shiftJacBit::setDist");
  if(ng_to < 0) return;

  rev >>= (ng_from - mp_to);

  for (int32_t step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    int32_t xy((fwd & rev).count());
    chr.nc[step] = xy/(xysum-xy);
  }
}

void shiftHamming::setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  rev <<= mp_from;
  for (int32_t step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    chr.mp[step] = (fwd ^ rev).count();
  }
  rev >>= (ng_from - mp_to);

  for (int32_t step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    chr.nc[step] = (fwd ^ rev).count();
  }
}

std::vector<int8_t> genVector(const std::vector<Read> &vReadref, const int32_t start, const int32_t end)
{
  std::vector<int8_t> array(end-start, 0);
  for (auto &x: vReadref) {
    if (!x.duplicate && my_range(x.F3, start, end-1)) ++array[x.F3 - start];
  }
  return array;
}

boost::dynamic_bitset<> genBitset(const std::vector<Read> &vReadref, const int32_t start, const int32_t end)
{
  boost::dynamic_bitset<> array(end-start);
  for (auto &x: vReadref) {
    if (!x.duplicate && my_range(x.F3, start, end-1))
      array.set(x.F3 - start);
  }
  return array;
}

namespace {
  void setSSPstats(SSPstats &p, const double bu, const double nsc, const double rlsc, const double rsc)
  {
    p.setnsc(nsc);
    p.setrlsc(rlsc);
    p.setrsc(rsc);
    p.setbu(bu);
  }
}


template <class T>
void genThread(T &dist, const SeqStatsGenomeSSP &genome, uint32_t chr_s, uint32_t chr_e, const std::string &prefix, const bool output_eachchr, const int32_t ng_to)
{
  for (uint32_t i=chr_s; i<=chr_e; ++i) {
    if (static_cast<int32_t>(genome.chr[i].getlen()) < ng_to) {
      std::cerr << "\nWarning: length of chromosome " << genome.chr[i].getname() << ": " << genome.chr[i].getlen()
		<< " is shorter than background distance " << ng_to << ". Skipped." << std::endl;
      //      std::cerr << "please specify shorter length with --ng_from and --ng_to options." << std::endl;
      continue;
    }
    std::cout << genome.chr[i].getname() << ".." << std::flush;

    dist.execchr(genome, i);
    dist.chr[i].setflen(dist.name);

    if (output_eachchr) {
      std::string filename = prefix + "." + genome.chr[i].getname() + ".csv";
      dist.outputmpChr(filename, i);
    }
  }
}

template <class T, class S>
void makeProfile(SSPstats &sspst, S &genome, const std::string &head, const std::string &typestr)
{
  DEBUGprint("makeProfile: " + typestr);
  T dist(sspst, genome);
  dist.printStartMessage();

  boost::thread_group agroup;

  std::string prefix(head + "." + typestr);
  if (typestr == "hdp" || typestr == "jaccard") {
    for (size_t i=0; i<genome.vsepchr.size(); ++i) {
      agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(genome), genome.vsepchr[i].s, genome.vsepchr[i].e, boost::cref(prefix), sspst.isEachchr(), sspst.getNgTo()));
    }
    agroup.join_all();
  } else {
    genThread(dist, genome, 0, genome.chr.size()-1, prefix, sspst.isEachchr(), sspst.getNgTo());
  }

  for (size_t i=0; i<genome.chr.size(); ++i) {
    if (genome.chr[i].isautosome()) dist.addmp2genome(i);
  }

  dist.setflen(dist.name);
  genome.dflen.setflen_ssp(dist.getnsci());

  if(sspst.getNgTo() < 0) return;

  std::string prefix2 = head + "." + typestr;
  dist.outputmpGenome(prefix2);

  if (typestr == "jaccard") setSSPstats(sspst, dist.getbackgroundUniformity(), dist.getnsc(), dist.getrlsc(), dist.getrsc());

  DEBUGprint("makeProfile: " + typestr + " done.");
  return;
}

void strShiftProfile(SSPstats &sspst, SeqStatsGenomeSSP &genome,
		     const std::string &head, const std::string &typestr)
{
  DEBUGprint("strShiftProfile...");

  if (typestr=="exjaccard")    makeProfile<shiftJacVec>(sspst, genome, head, typestr);
  else if (typestr=="jaccard") makeProfile<shiftJacBit>(sspst, genome, head, typestr);
  else if (typestr=="ccp")     makeProfile<shiftCcp>(sspst, genome, head, typestr);
  else if (typestr=="hdp")     makeProfile<shiftHamming>(sspst, genome, head, typestr);
/*  else if (typestr=="drompa") {
    if (genome.dflen.isallchr()) makeProfile<shiftJacBit>(sspst, genome, head, "jaccard");
    else makeProfile_forDROMPA<shiftJacBit>(sspst, genome, head, "jaccard");
  }*/

  DEBUGprint("strShiftProfile done.");
  return;
}
