/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#include "ssp_shiftprofile.h"
#include "ssp_shiftprofile_p.h"
#include "macro.h"
#include <map>
#include <boost/thread.hpp>

void addmp(std::map<int32_t, double> &mpto, const std::map<int32_t, double> &mpfrom, double w)
{
  for(auto itr = mpfrom.begin(); itr != mpfrom.end(); ++itr) {
    if(!std::isnan(itr->second)) mpto[itr->first] += itr->second * w;
  }
}

double getJaccard(int32_t step, int32_t to, int32_t xysum, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  int32_t xy(0);
  for(int32_t j=mp_from; j<to; ++j) if(fwd[j] * rev[j+step]) xy += std::max(fwd[j], rev[j+step]);
  return (xy/static_cast<double>(xysum-xy));
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
  moment<int8_t> x(fwd, mp_from, chr.width-ng_to);
  moment<int8_t> y(rev, mp_from, chr.width-ng_to);

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

  double val = 1/(x.getsd() * y.getsd() * (chr.width-ng_to - mp_from - 1));
  for(auto itr = chr.mp.begin(); itr != chr.mp.end(); ++itr) itr->second *= val;
  for(auto itr = chr.nc.begin(); itr != chr.nc.end(); ++itr) itr->second *= val;
}

void shiftJacBit::setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  int32_t xysum(fwd.count() + rev.count());

  rev <<= mp_from;
  for(int32_t step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    int32_t xy((fwd & rev).count());
    //    chr.mp[step] = xy;
    chr.mp[step] = xy/static_cast<double>(xysum-xy);
  }
  rev >>= (ng_from - mp_to);

  for(int32_t step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    int32_t xy((fwd & rev).count());
    chr.nc[step] = xy/static_cast<double>(xysum-xy);
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

std::vector<int8_t> genVector(const strandData &seq, int32_t start, int32_t end)
{
  std::vector<int8_t> array(end-start, 0);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1)) ++array[x.F3 - start];
  }
  return array;
}

std::vector<int8_t> shiftFragVar::genVector4FixedReadsNum(const strandData &seq, int32_t start, int32_t end)
{
  std::vector<int8_t> array(end-start, 0);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1)){
      if(rand() >= r4cmp) continue;
      ++array[x.F3 - start];
      ++numUsed4FCS;
    }
  }
  return array;
}

boost::dynamic_bitset<> genBitset(const strandData &seq, int32_t start, int32_t end)
{
  boost::dynamic_bitset<> array(end-start);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1))
      array.set(x.F3 - start);
  }
  return array;
}

template <class T>
void genThread(T &dist, const Mapfile &p, uint32_t chr_s, uint32_t chr_e, std::string typestr, const bool output_eachchr) {
  for(uint32_t i=chr_s; i<=chr_e; ++i) {
    std::cout << p.genome.chr[i].name << ".." << std::flush;

    dist.execchr(p, i);
    dist.chr[i].setflen(dist.name);
    
    std::string filename = p.getprefix() + "." + typestr + "." + p.genome.chr[i].name + ".csv";
    if(output_eachchr) dist.outputmpChr(filename, i);
  }
}

template <class T>
void makeProfile(Mapfile &p, const std::string &typestr, const MyOpt::Variables &values)
{
  T dist(p, values);
  dist.printStartMessage();
  
  boost::thread_group agroup;
  boost::mutex mtx;

  if(typestr == "hdp" || typestr == "jaccard") {
    //agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), 0, 0, typestr, values.count("eachchr")));
    for(uint32_t i=0; i<p.genome.vsepchr.size(); i++) {
      agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), p.genome.vsepchr[i].s, p.genome.vsepchr[i].e, typestr, values.count("eachchr")));
    }
    agroup.join_all();
  } else {
    genThread(dist, p, 0, p.genome.chr.size()-1, typestr, values.count("eachchr"));
    //genThread(dist, p, 0, 0, typestr);
  }

  // set fragment length;
  for(uint32_t i=0; i<p.genome.chr.size(); ++i) {
    if(p.genome.chr[i].isautosome()) dist.addmp2genome(i);
  }

  dist.setflen(dist.name);
  p.seteflen(dist.getnsci());

  std::string prefix = p.getprefix() + "." + typestr;
  dist.outputmpGenome(prefix);

  if(typestr == "jaccard") {
    p.setbackgroundUniformity(dist.getbackgroundUniformity());
  }

  return;
}

void strShiftProfile(const MyOpt::Variables &values, Mapfile &p, std::string typestr)
{
  if(typestr=="exjaccard")    makeProfile<shiftJacVec>(p, typestr, values);
  else if(typestr=="jaccard") makeProfile<shiftJacBit>(p, typestr, values);
  else if(typestr=="ccp")     makeProfile<shiftCcp>(p, typestr, values);
  else if(typestr=="hdp")     makeProfile<shiftHamming>(p, typestr, values);
  
  return;
}

void makeRscript(const std::string prefix)
{
  std::string Rscript(prefix + ".FCS.R");
  std::ofstream out(Rscript);
  out << "data <- read.csv('" << prefix << ".acfp.csv', header=TRUE, row.names=1, sep='\t', quote='')" << std::endl;
  out << "colname <- colnames(data)" << std::endl;
  out << "colnames(data) <- colname[-1]" << std::endl;
  out << "data <- data[,-11]" << std::endl;
  out << "pdf('" << prefix << ".FCS.pdf', height=7, width=14)" << std::endl;
  out << "par(mfrow=c(1,2))" << std::endl;
  out << "cols <- rainbow(10)" << std::endl;
  out << "plot(0, 0, type = 'n', xlim = range(1:nrow(data)), ylim = range(data), xlab = 'Neighboring distance (bp)', ylab = 'Accumulated prop')" << std::endl;
  out << "for (i in 1:ncol(data)) { lines(1:nrow(data), data[,i], col=cols[i])}" << std::endl;
  out << "legend('bottomright', legend = colnames(data), lty = 1, col = cols)" << std::endl;
  out << "data <- read.csv('" << prefix << ".fcs.csv', header=TRUE, skip=2, sep='\t', quote='')" << std::endl;
  out << "plot(data[,1],data[,2], log='x', type='l', xlab = 'Read-pair distance (bp)', ylab = 'Fragment cluster score')" << std::endl;
  out << "dev.off()" << std::endl;

  std::string command = "R --vanilla < " + Rscript + " > " + Rscript + ".log 2>&1";
  std::cout << command << std::endl;
  int32_t return_code = system(command.c_str());
  if(WEXITSTATUS(return_code)) {
    std::cerr << "Warning: command " << command << "return nonzero status." << std::endl;
  }
  
  return;
}

void makeFCSProfile(const MyOpt::Variables &values, Mapfile &p, const std::string &typestr)
{
  shiftFragVar dist(p, values, p.getflen(values));
  dist.printStartMessage();

  for(uint32_t i=0; i<p.genome.chr.size(); ++i) {
    if(p.genome.chr[i].isautosome()) {
      std::cout << p.genome.chr[i].name << ".." << std::flush;
      dist.execchr(p, i);
      if(values.count("eachchr")) {
	std::string filename = p.getprefix() + "." + typestr + "." + p.genome.chr[i].name + ".csv";
	dist.outputmpChr(filename, i);
      }
    }
  }
  std::cout << "\nread number for calculating FCS: " << dist.getnumUsed4FCS() << std::endl;

  std::string filename1 = p.getprefix() + ".acfp.csv";
  dist.printacfp(filename1);
  std::string filename2 = p.getprefix() + "." + typestr + ".csv";
  dist.outputmpGenome(filename2);

  makeRscript(p.getprefix());

  return;
}

void genThreadFragVar(ReadShiftProfile &chr, std::map<int32_t, FragmentVariability> &acfp, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev, const std::vector<double> &fvback, const int32_t s, const int32_t e, boost::mutex &mtx)
{
  for(int32_t step=s; step<e; ++step) {
    FragmentVariability fv;
    fv.setVariability(step, chr.start, chr.end, fwd, rev);

    double diffMax(0);
    for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) {
      //      std::cout << fv.getAccuOfDistanceOfFragment(k) << "\t" << fvback[k] << std::endl;
      diffMax = std::max(diffMax, fv.getAccuOfDistanceOfFragment(k) - fvback[k]);
    }
    chr.setmp(step, diffMax, mtx);
    acfp[step].add2genome(fv, mtx);
  }
}

void shiftFragVar::setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  boost::thread_group agroup;
  boost::mutex mtx;

  std::vector<double> fvback(sizeOfvDistOfDistaneOfFrag,0);
  int32_t n(0);
  for(int32_t step=ng_from_fcs; step<ng_to_fcs; step+=ng_step_fcs) {
    FragmentVariability fv;
    fv.setVariability(step, chr.start, chr.end, fwd, rev);    
    for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) {
      fvback[k] += fv.getAccuOfDistanceOfFragment(k);
    }
    ++n;
  }
  for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) fvback[k] /= n;

  std::vector<int32_t> v{flen, chr.getlenF3()};
  std::copy(v4acfp.begin(), v4acfp.end(), std::back_inserter(v));
  for(auto x: v) {
    FragmentVariability fv;
    fv.setVariability(x, chr.start, chr.end, fwd, rev);
    
    double diffMax(0);
    for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) {
      //      std::cout << fv.getAccuOfDistanceOfFragment(k) << "\t" << fvback[k] << std::endl;
      diffMax = std::max(diffMax, fv.getAccuOfDistanceOfFragment(k) - fvback[k]);
    }
    chr.setmp(x, diffMax, mtx);
    acfp[x].add2genome(fv, mtx);
  }
  
  return;
}
