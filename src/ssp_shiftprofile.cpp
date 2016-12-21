/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#include "ssp_shiftprofile.h"
#include "ssp_shiftprofile_p.h"
#include "macro.h"
#include <map>
#include <boost/thread.hpp>

void addmp(std::map<int, double> &mpto, const std::map<int, double> &mpfrom, double w)
{
  for(auto itr = mpfrom.begin(); itr != mpfrom.end(); ++itr) {
    mpto[itr->first] += itr->second * w;
  }
}

double getJaccard(int step, int width, int xysum, const std::vector<char> &fwd, const std::vector<char> &rev)
{
  int xy(0);
  for(int j=mp_from; j<width-ng_to; ++j) if(fwd[j] * rev[j+step]) xy += std::max(fwd[j], rev[j+step]);
  return (xy/static_cast<double>(xysum-xy));
}

void genThreadJacVec(ReadShiftProfile &chr, int xysum, const std::vector<char> &fwd, const std::vector<char> &rev, int s, int e, boost::mutex &mtx)
{
  for(int step=s; step<e; ++step) {
    chr.setmp(step, getJaccard(step, chr.width, xysum, fwd, rev), mtx);
  }
}

void shiftJacVec::setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev)
{
  int xx = accumulate(fwd.begin(), fwd.end(), 0);
  int yy = accumulate(rev.begin(), rev.end(), 0);

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(&genThreadJacVec, boost::ref(chr), xx+yy, boost::cref(fwd), boost::cref(rev), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for(int step=ng_from; step<ng_to; step+=ng_step) {
    chr.nc[step] = getJaccard(step, chr.width-ng_to, xx+yy, fwd, rev);
  }
}

void genThreadCcp(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev, double mx, double my, const int s, const int e, boost::mutex &mtx)
{
  for(int step=s; step<e; ++step) {
    double xy(0);
    for(int j=mp_from; j<chr.width-ng_to; ++j) xy += (fwd[j] - mx) * (rev[j+step] - my);
    chr.setmp(step, xy, mtx);
  }
}

void shiftCcp::setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev)
{
  moment<char> x(fwd, mp_from, chr.width-ng_to);
  moment<char> y(rev, mp_from, chr.width-ng_to);

  boost::thread_group agroup;
  boost::mutex mtx;
  for(uint i=0; i<seprange.size(); i++) {
    agroup.create_thread(bind(genThreadCcp, boost::ref(chr), boost::cref(fwd), boost::cref(rev), x.getmean(), y.getmean(), seprange[i].start, seprange[i].end, boost::ref(mtx)));
  }
  agroup.join_all();

  for(int step=ng_from; step<ng_to; step+=ng_step) {
    double xy(0);
    for(int j=mp_from; j<chr.width-ng_to; ++j) xy += (fwd[j] - x.getmean()) * (rev[j+step] - y.getmean());
    chr.nc[step] = xy;
  }

  double val = 1/(x.getsd() * y.getsd() * (chr.width-ng_to - mp_from - 1));
  for(auto itr = chr.mp.begin(); itr != chr.mp.end(); ++itr) itr->second *= val;
  for(auto itr = chr.nc.begin(); itr != chr.nc.end(); ++itr) itr->second *= val;
}

void shiftJacBit::setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  int xysum(fwd.count() + rev.count());

  rev <<= mp_from;
  for(int step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    int xy((fwd & rev).count());
    //    chr.mp[step] = xy;
    chr.mp[step] = xy/static_cast<double>(xysum-xy);
  }
  rev >>= (ng_from - mp_to);

  for(int step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    int xy((fwd & rev).count());
    chr.nc[step] = xy/static_cast<double>(xysum-xy);
  }
}

void shiftHamming::setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev)
{
  rev <<= mp_from;
  for(int step=-mp_from; step<mp_to; ++step) {
    rev >>= 1;
    chr.mp[step] = (fwd ^ rev).count();
  }
  rev >>= (ng_from - mp_to);
  
  for(int step=ng_from; step<ng_to; step+=ng_step) {
    rev >>= ng_step;
    chr.nc[step] = (fwd ^ rev).count();
  }
}

std::vector<char> genVector(const strandData &seq, int start, int end)
{
  std::vector<char> array(end-start, 0);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1))
      ++array[x.F3 - start];
  }
  return array;
}

std::vector<char> genVector4FixedReadsNum(const strandData &seq, int start, int end, const double r4cmp)
{
  std::vector<char> array(end-start, 0);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1)){
      if(rand() >= r4cmp) continue;
      ++array[x.F3 - start];
    }
  }

  return array;
}

boost::dynamic_bitset<> genBitset(const strandData &seq, int start, int end)
{
  boost::dynamic_bitset<> array(end-start);
  for (auto x: seq.vRead) {
    if(!x.duplicate && RANGE(x.F3, start, end-1))
      array.set(x.F3 - start);
  }
  return array;
}

template <class T>
void genThread(T &dist, const Mapfile &p, uint chr_s, uint chr_e, std::string typestr, const bool output_eachchr) {
  for(uint i=chr_s; i<=chr_e; ++i) {
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
  T dist(p, values["threads"].as<int>());
  dist.printStartMessage();
  
  boost::thread_group agroup;
  boost::mutex mtx;

  if(typestr == "hdp" || typestr == "jaccard") {
    //agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), 0, 0, typestr, values.count("output_eachchr")));
    for(uint i=0; i<p.genome.vsepchr.size(); i++) {
      agroup.create_thread(bind(genThread<T>, boost::ref(dist), boost::cref(p), p.genome.vsepchr[i].s, p.genome.vsepchr[i].e, typestr, values.count("output_eachchr")));
    }
    agroup.join_all();
  } else {
    genThread(dist, p, 0, p.genome.chr.size()-1, typestr, values.count("output_eachchr"));
    //genThread(dist, p, 0, 0, typestr);
  }

  // set fragment length;
  for(uint i=0; i<p.genome.chr.size(); ++i) {
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

  std::string command = "R --vanilla < " + Rscript + " | tee " + Rscript + ".log";
  int return_code = system(command.c_str());
  
  return;
}

void makeFCSProfile(const MyOpt::Variables &values, Mapfile &p, const std::string &typestr)
{
  shiftFragVar dist(p, values["threads"].as<int>(), p.getflen(values), values.count("fcsfull"));
  dist.printStartMessage();

  int numRead4fcs(values["nfcs"].as<int>());
  double r(1);
  //  std::cout << "backgrounduniformity: " << p.getbackgroundUniformity() << std::endl;
  if(numRead4fcs) r = numRead4fcs/static_cast<double>(dist.getnread());
  
  if(r>1){
    std::cerr << "\nWarning: number of reads for Fragment variability is < "<< (int)(numRead4fcs/NUM_1M) <<" million.\n";
    dist.lackOfReads_on();
  }

  double r4cmp = r*RAND_MAX;

  for(uint i=0; i<=p.genome.chr.size()-1; ++i) {
    if(p.genome.chr[i].isautosome()) {
      std::cout << p.genome.chr[i].name << ".." << std::flush;
      dist.execchr(p, i, r4cmp);
      std::string filename = p.getprefix() + "." + typestr + "." + p.genome.chr[i].name + ".csv";
      if(values.count("output_eachchr")) dist.outputmpChr(filename, i);
      dist.addmp2genome(i);
    }
  }

  std::string filename1 = p.getprefix() + ".acfp.csv";
  dist.printacfp(filename1);
  std::string filename2 = p.getprefix() + "." + typestr + ".csv";
  dist.outputmpGenome(filename2);

  makeRscript(p.getprefix());

  return;
}

void genThreadFragVar(ReadShiftProfile &chr, std::map<int, FragmentVariability> &acfp, const std::vector<char> &fwd, const std::vector<char> &rev, const std::vector<double> &fvback, const int s, const int e, boost::mutex &mtx)
{
  for(int step=s; step<e; ++step) {
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

void shiftFragVar::setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev)
{
  boost::thread_group agroup;
  boost::mutex mtx;

  std::vector<double> fvback(sizeOfvDistOfDistaneOfFrag,0);
  int n(0);
  //  for(int step=ng_from; step<ng_to; step+=ng_step) {
  for(int step=NUM_100K; step<NUM_1M; step+=NUM_100K) {
    FragmentVariability fv;
    fv.setVariability(step, chr.start, chr.end, fwd, rev);    
    for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) {
      fvback[k] += fv.getAccuOfDistanceOfFragment(k);
    }
    ++n;

  }
  for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag; ++k) fvback[k] /= n;

  if (fcsfull) {
    for(uint i=0; i<seprange.size(); i++) {
      agroup.create_thread(bind(&genThreadFragVar, boost::ref(chr), boost::ref(acfp), boost::cref(fwd), boost::cref(rev), boost::cref(fvback), seprange[i].start, seprange[i].end, boost::ref(mtx)));
    }
    agroup.join_all();
  } else {
    std::vector<int> v{flen, chr.getlenF3()};
    std::copy(v4acfp.begin(), v4acfp.end(), std::back_inserter(v));
    for(auto x: v) {
    //    for(auto x: v4acfp) {
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
  }

  return;
}

/*void scanRepeatRegion(const std::vector<char> &fwd, const std::vector<char> &rev)
{
  int SizeOfFragOverlapDist(10000);
  std::vector<int> FragOverlapDist(SizeOfFragOverlapDist,0);

  int size(fwd.size());
  std::vector<short> array(size, 0);
  int fraglen=1000;
  int last(0);
  for(int i=0; i<size - fraglen; ++i) {
    if(fwd[i] && rev[i+fraglen]) {
      for(int j=0; j<fraglen; ++j) ++array[i+j];
      if(i-last==1) std::cout << last << "-" << i << std::endl;
      last=i;
    }
  }

  //  for(int i=0; i<size; ++i) if(array[i]>1) std::cout << i << "\t" << array[i] << std::endl;
  return;
}
*/
