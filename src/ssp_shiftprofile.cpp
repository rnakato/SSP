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
  MyStatistics::moment<int8_t> x(fwd, mp_from, chr.width-ng_to);
  MyStatistics::moment<int8_t> y(rev, mp_from, chr.width-ng_to);

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
    chr.mp[step] = getratio(xy, xysum-xy);
  }
  rev >>= (ng_from - mp_to);

  for(int32_t step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    int32_t xy((fwd & rev).count());
    chr.nc[step] = getratio(xy, xysum-xy);
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

std::vector<int8_t> genVector(const std::vector<Read> &vReadref, int32_t start, int32_t end)
{
  std::vector<int8_t> array(end-start, 0);
  for (auto &x: vReadref) {
    if(!x.duplicate && my_range(x.F3, start, end-1)) ++array[x.F3 - start];
  }
  return array;
}

std::vector<int8_t> shiftFragVar::genVector4FixedReadsNum(const std::vector<Read> &vReadref, int32_t start, int32_t end) {
  std::vector<int8_t> array(end-start, 0);
  for (auto &x: vReadref) {
    if(!x.duplicate && my_range(x.F3, start, end-1)){
      if(rand() >= r4cmp) continue;
      ++array[x.F3 - start];
      ++numUsed4FCS;
    }
  }
  return array;
}

boost::dynamic_bitset<> genBitset(const std::vector<Read> &vReadref, int32_t start, int32_t end)
{
  boost::dynamic_bitset<> array(end-start);
  for (auto &x: vReadref) {
    if(!x.duplicate && my_range(x.F3, start, end-1))
      array.set(x.F3 - start);
  }
  return array;
}

namespace {
  void setSSPstats(SSPstats &p, const double bu, const double nsc, const double rsc) {
    p.setnsc(nsc);
    p.setrsc(rsc);
    p.setbu(bu);
  }
  
  void setFCSstats(SSPstats &p, const double fcsread, const double fcsflen, const double fcs1k, const double fcs10k, const double fcs100k) {
    p.setfcsread(fcsread);
    p.setfcsflen(fcsflen);
    p.setfcs1k(fcs1k);
    p.setfcs10k(fcs10k);
    p.setfcs100k(fcs100k);
  }
}

template <class T>
void genThread(T &dist, const SeqStatsGenome &genome, uint32_t chr_s, uint32_t chr_e, const std::string &prefix, const bool output_eachchr) {
  for(uint32_t i=chr_s; i<=chr_e; ++i) {
    std::cout << genome.chr[i].getname() << ".." << std::flush;

    dist.execchr(genome, i);
    dist.chr[i].setflen(dist.name);
    
    std::string filename = prefix + "." + genome.chr[i].getname() + ".csv";
    if(output_eachchr) dist.outputmpChr(filename, i);
  }
}

template <class T>
void makeProfile(SSPstats &sspst, SeqStatsGenome &genome, const std::string &head, const std::string &typestr, const MyOpt::Variables &values)
{
  T dist(genome, values);
  dist.printStartMessage();
  
  boost::thread_group agroup;
  boost::mutex mtx;

  std::string prefix(head + "." + typestr);
  if(typestr == "hdp" || typestr == "jaccard") {
    for(size_t i=0; i<genome.vsepchr.size(); i++) {
      agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(genome), genome.vsepchr[i].s, genome.vsepchr[i].e, boost::cref(prefix), values.count("eachchr")));
    }
    agroup.join_all();
  } else {
    genThread(dist, genome, 0, genome.chr.size()-1, prefix, values.count("eachchr"));
  }

  for(size_t i=0; i<genome.chr.size(); ++i) {
    if(genome.chr[i].isautosome()) dist.addmp2genome(i);
  }

  dist.setflen(dist.name);
  genome.dflen.setflen_ssp(dist.getnsci());

  std::string prefix2 = head + "." + typestr;
  dist.outputmpGenome(prefix2);

  if(typestr == "jaccard") setSSPstats(sspst, dist.getbackgroundUniformity(), dist.getnsc(), dist.getrsc());

  return;
}

void strShiftProfile(SSPstats &sspst, const MyOpt::Variables &values, SeqStatsGenome &genome, const std::string &head, const std::string &typestr)
{
  if(typestr=="exjaccard")    makeProfile<shiftJacVec>(sspst, genome, head, typestr, values);
  else if(typestr=="jaccard") makeProfile<shiftJacBit>(sspst, genome, head, typestr, values);
  else if(typestr=="ccp")     makeProfile<shiftCcp>(sspst, genome, head, typestr, values);
  else if(typestr=="hdp")     makeProfile<shiftHamming>(sspst, genome, head, typestr, values);

  return;
}

void makeRscript(const std::string prefix)
{
  std::string Rscript(prefix + ".FCS.R");
  std::ofstream out(Rscript);
  out << "data <- read.csv('" << prefix << ".pnf.csv', header=TRUE, row.names=1, sep='\t', quote='')" << std::endl;
  out << "colnames(data) <- colnames(data)[-1]" << std::endl;
  out << "data <- data[,-ncol(data)]" << std::endl;
  out << "ncol <- ncol(data)/2" << std::endl;
  out << "nrow <- nrow(data)" << std::endl;
  out << "cols <- rainbow(ncol)" << std::endl;
  out << "cpnf <- data[,(ncol+1):(ncol*2)]" << std::endl;
  out << "x <- seq(1, nrow, 10)" << std::endl;
  out << "cpnf10 <- cpnf[x,]" << std::endl;
  out << "pnf <- rbind(cpnf10[-1,],cpnf[nrow,]) - cpnf10" << std::endl;
  out << "pdf('" << prefix << ".FCS.pdf', height=5, width=15)" << std::endl;
  out << "par(mfrow=c(1,3))" << std::endl;
  // Proportion of NN fragments
  out << "plot(0, 0, type = 'n', xlim = range(1:nrow), ylim = range(pnf), xlab = 'Neighboring distance (bp)', ylab = 'Proportion of nearest neibor fragments')" << std::endl;
  out << "for (i in 1:ncol) { lines(x, pnf[,i], col=cols[i])}" << std::endl;
  out << "legend('bottomright', legend = colnames(pnf), lty = 1, col = cols)" << std::endl;
  // Cumurative proportion
  out << "plot(0, 0, type = 'n', xlim = range(1:nrow), ylim = range(cpnf), xlab = 'Neighboring distance (bp)', ylab = 'Cumulative proportion')" << std::endl;
  out << "for (i in 1:ncol) { lines(1:nrow, cpnf[,i], col=cols[i])}" << std::endl;
  out << "legend('bottomright', legend = colnames(cpnf), lty = 1, col = cols)" << std::endl;
  // FCS
  out << "data <- read.csv('" << prefix << ".fcs.csv', header=TRUE, skip=4, sep='\t', quote='')" << std::endl;
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

void makeFCSProfile(SSPstats &sspst, const MyOpt::Variables &values, const SeqStatsGenome &genome, const std::string &head, const std::string &typestr)
{
  shiftFragVar dist(genome, values, genome.dflen.getflen());
  dist.printStartMessage();

  for(uint32_t i=0; i<genome.chr.size(); ++i) {
    if(!genome.chr[i].isautosome()) continue;
    std::cout << genome.chr[i].getname() << ".." << std::flush;
    dist.execchr(genome, i);
    if(values.count("eachchr")) {
      std::string filename = head + "." + typestr + "." + genome.chr[i].getname() + ".csv";
      dist.outputfcsChr(filename, i);
    }
  }
  std::cout << "\nread number for calculating FCS: " << dist.getnumUsed4FCS() << std::endl;

  std::string filename1 = head + ".pnf.csv";
  dist.printdistpnf(filename1);
  std::string filename2 = head + "." + typestr + ".csv";
  dist.outputfcsGenome(filename2);

  makeRscript(head);

  setFCSstats(sspst, dist.getMPread(), dist.getMPflen(), dist.getMP1k(), dist.getMP10k(), dist.getMP100k());

  return;
}

double getFCS(PropNeighborFrag &fv, std::vector<double> &fvbg)
{
  double diffMax(0);
  for(size_t k=0; k<sizeOfvNeighborFrag; ++k) {
    diffMax = std::max(diffMax, fv.getCumulativePNF(k) - fvbg[k]);
  }
  return diffMax;
}

void shiftFragVar::setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev)
{
  boost::thread_group agroup;
  boost::mutex mtx;

  // make fv for background
  std::vector<double> fvbg(sizeOfvNeighborFrag, 0);
  int32_t n(0);
  for(int32_t flen=ng_from_fcs; flen<ng_to_fcs; flen += ng_step_fcs) {
    PropNeighborFrag fv;
    fv.setNeighborFrag(flen, chr.start, chr.end, fwd, rev);    
    for(size_t k=0; k<sizeOfvNeighborFrag; ++k) {
      fvbg[k] += fv.getCumulativePNF(k);
    }
    ++n;
  }
  for(size_t k=0; k<sizeOfvNeighborFrag; ++k) fvbg[k] /= n;

  // calculate FCS for each step
  std::vector<int32_t> v{flen, chr.getlenF3()};
  std::copy(v4pnf.begin(), v4pnf.end(), std::back_inserter(v));
  for(auto len: v) {
    PropNeighborFrag fv;
    fv.setNeighborFrag(len, chr.start, chr.end, fwd, rev);

    double fcs = getFCS(fv, fvbg);
    chr.setmp(len, fcs, mtx);
    distpnf[len].add2genome(fv, mtx);
  }
  return;
}
