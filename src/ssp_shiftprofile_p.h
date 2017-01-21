/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#ifndef _SSP_SHIFTPROFILE_P_H_
#define _SSP_SHIFTPROFILE_P_H_

#include "pw_gv.h"
#include <boost/dynamic_bitset.hpp>
#include "alglib/alglib.h"

namespace {
  const int32_t mp_from(500);
  const int32_t mp_to(1500);
  const int32_t sizeOfvNeighborFrag = 5000;
  const std::vector<int32_t> v4pnf{50, 100, 150, 500, 1000, 2000, 3000, 10000, 100000, 1000000};
}

std::vector<int8_t> genVector(const strandData &seq, int32_t start, int32_t end);
std::vector<int8_t> genVector4FixedReadsNum(const strandData &seq, int32_t start, int32_t end, const double r4cmp, uint32_t &);
boost::dynamic_bitset<> genBitset(const strandData &seq, int32_t, int32_t);
void addmp(std::map<int32_t, double> &, const std::map<int32_t, double> &, double w);

double getmpmean(const std::map<int32_t, double> mp, int32_t s, int32_t e) {
  double m(0);
  int32_t n(0);
  if(s > e) {
    std::cerr << "Error: s " << s << "> e " << e <<"for getmpmean" << std::endl;
    return -1;
  }
  for(int32_t i=s; i <= e; ++i) {
    m += mp.at(i);
    ++n;
  }
  return m/n;
}

class PropNeighborFrag {
  double sumOfvNeighborFrag;
  std::vector<int32_t> vNeighborFrag;

 public:
 PropNeighborFrag(): sumOfvNeighborFrag(0), vNeighborFrag(sizeOfvNeighborFrag, 0) {}

  void setNeighborFrag(const int32_t flen, const int32_t start, const int32_t end,
		      const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev) {
    if(start < 0 || flen < 0) {
      std::cerr << "error: invalid start " << start << "or flen " << flen << "for setNeighborFrag." << std::endl;
      return;
    }
    
    int32_t last(start);
    for(int32_t i=start; i<end-flen; ++i) {
      if(fwd[i] && rev[i+flen]) {
	int32_t distance = i-last;
	if(distance <= 0) {
	  //	  std::cerr << "error: invalid distance " << distance << std::endl;
	}else if(distance < sizeOfvNeighborFrag-1) ++vNeighborFrag[distance];
	else ++vNeighborFrag[sizeOfvNeighborFrag-1];
	last = i;
      }
    }
    sumOfvNeighborFrag = accumulate(vNeighborFrag.begin(), vNeighborFrag.end(), 0);
  }
  
  double getPNF(const int32_t i) const {
    if(i<0 || i>= sizeOfvNeighborFrag) {
      std::cerr << "error: invalid num " << i << "for getPNF max: " << sizeOfvNeighborFrag << std::endl;
      return -1;
    }
    else {
      return sumOfvNeighborFrag ? vNeighborFrag[i] / sumOfvNeighborFrag : 0;
    }
  }
  double getCumulativePNF(const int32_t i) const {
    if(i<0 || i>= sizeOfvNeighborFrag) {
      std::cerr << "error: invalid num " << i << "for getCumulativePNF max: " << sizeOfvNeighborFrag << std::endl;
      return -1;
    }
    else {
      double cpnf(0);
      for(int32_t j=0; j<=i; ++j) cpnf += getPNF(j);
      return cpnf;
    }
  }
  void add2genome(const PropNeighborFrag &x, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    for(size_t i=0; i<vNeighborFrag.size(); ++i) vNeighborFrag[i] += x.vNeighborFrag[i];
    sumOfvNeighborFrag += x.sumOfvNeighborFrag;
  }
};

class ReadShiftProfile {
  int32_t lenF3;
  double r;
  double bk;
  int32_t bk_from;

 protected:
  double nsc;
  double rsc;
  double rlsc;
  int32_t nsci;
  uint64_t len;
  uint64_t nread;
  uint32_t num4ssp;
  double backgroundUniformity;
  
 public:
  std::map<int32_t, double> mp;
  std::map<int32_t, double> nc;
  int32_t start;
  int32_t end;
  int32_t width;

  double rchr;

 ReadShiftProfile(const int32_t lenf3, const int32_t b, const int32_t n4s, int32_t s=0, int32_t e=0, int64_t n=0, int64_t l=0):
  lenF3(lenf3), r(0), bk(0), bk_from(b), nsc(0), rsc(0), rlsc(0), nsci(0), len(l), nread(n), num4ssp(n4s), backgroundUniformity(0), start(s), end(e), width(e-s), rchr(1) {}
  virtual ~ReadShiftProfile() {}
  void setmp(const int32_t i, const double val, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    mp[i] = val;
  }

  double getbackgroundUniformity() const { return backgroundUniformity; }
  void setrchr(const uint64_t n) { rchr = n ? nread/static_cast<double>(n): 0; }
  int32_t getlenF3() const { return lenF3; }
  double getnsc() const { return nsc; }
  double getrsc() const { return rsc; }
  int32_t getnsci() const { return nsci; }
  uint64_t getnread() const { return nread; }
  uint64_t getlen() const { return len; }
  double getmpsum() const {
    double sum(0);
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) sum += itr->second;
    return sum;
  }
  void setControlRatio() {
    int32_t n(0);
    for(auto itr = nc.begin(); itr != nc.end(); ++itr) {
      if(itr->first >= bk_from) {
	bk += itr->second;
	++n;
      }
    }
    bk /= n;
    r = 1/bk;
  }

  void setflen(const std::string &name) {
    int32_t threwidth(5);
    int32_t leftend(lenF3*1.2);
    if(leftend>150) leftend=150;
    //    if(leftend<lenF3) leftend=lenF3;

    setControlRatio();
    nsc = mp[mp_to-1];

    std::map<int32_t, double> mpsmooth;
    for(int32_t i=mp_to-1-2; i > leftend-threwidth; --i) {
      mpsmooth[i] = getmpmean(mp, i-2, i+2);
    }

    for(int32_t i=mp_to-1-threwidth-2; i > leftend; --i) {
      if(name == "Hamming distance") {
	if(mpsmooth.at(i) < mpsmooth.at(i+threwidth) || mpsmooth.at(i) < mpsmooth.at(i-threwidth)) continue;
	if(nsc > mp.at(i)) {
	  nsc  = mp.at(i);
	  nsci = i;
	}
      } else {
	if(mpsmooth.at(i) < mpsmooth.at(i+threwidth) || mpsmooth.at(i) < mpsmooth.at(i-threwidth)) continue;
	double s(mp.at(i)*r);
	if(nsc < s) {
	  nsc  = s;
	  rsc  = (mp.at(i) - bk)/(mp.at(lenF3) - bk);
	  nsci = i;
	}
      }
    }
  }

  double getMPread() const { return  mp.at(lenF3); }
  double getMP1k()   const { return  mp.at(1000); }
  double getMP10k()  const { return  mp.at(10000); }
  double getMP100k()  const { return  mp.at(100000); }

  void print2file(const std::string &filename, const std::string &name) {
    if(!nread) {
      std::cerr << filename << ": no read" << std::endl;
    }
    double sum(getmpsum());
    double rRPKM = (num4ssp/static_cast<double>(nread)) / (NUM_100M/static_cast<double>(len));

    double be(bk * rRPKM);
    double const_bu = num4ssp/static_cast<double>(4*NUM_100M - num4ssp);  // 1/39 N/(4*L-N), N=10M, L=100M
    //    std::cout << "####### " << num4ssp << "\t  " << const_bu << "\t" << rRPKM << "\t" << bk << std::endl;
    //    std::cout << "####### " << nread << "\t  " << NUM_100M << "\t" << len << std::endl;
    rlsc = getMPread() *r;
    backgroundUniformity = const_bu / be;

    std::ofstream out(filename);
    out << "NSC\t" << nsc << std::endl;
    out << "RSC\t" << rsc << std::endl;
    out << "RLSC\t"<< rlsc << std::endl;
    out << "Estimated fragment length\t" << nsci << std::endl;
    out << "Background enrichment\t" << be << std::endl;
    out << "Background uniformity\t" << backgroundUniformity << std::endl;

    out << "Strand shift\t" << name << "\tProportion\tper " << num4ssp/NUM_1M << "M reads for 100Mbp len\tper control" << std::endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) 
      out << itr->first            << "\t"
	  << itr->second           << "\t"
	  << (itr->second/sum)     << "\t"
	  << (itr->second * rRPKM) << "\t"
	  << (itr->second * r)     << std::endl;
    for(auto itr = nc.begin(); itr != nc.end(); ++itr) 
      out << itr->first            << "\t"
	  << itr->second           << "\t"
	  << (itr->second/sum)     << "\t"
	  << (itr->second * rRPKM) << "\t"
	  << (itr->second * r)     << std::endl;
  }
  void print2file4fcs(const std::string &filename, const std::string &name, const int32_t flen, const bool lackOfReads) const {
    if(!nread) std::cerr << filename << ": no read" << std::endl;

    //    std::cout << flen <<"," << lenF3 <<std::endl;

    std::ofstream out(filename);
    std::string str("");
    if(lackOfReads) str = " (read number is insufficient)";
    out << "Read length"     << str << "\t" << getMPread() << std::endl;
    out << "Fragment length" << str << "\t" << mp.at(flen) << std::endl;
    out << "Broad (1 kbp)"   << str << "\t" << getMP1k()   << std::endl;
    out << "Broad (10 kbp)"  << str << "\t" << getMP10k()  << std::endl;
    out << "Strand shift\t"  << name << std::endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr)
      //     if(itr->first != flen && itr->first != lenF3)
	out << itr->first << "\t" << itr->second << std::endl;
  }
};

class ReadShiftProfileGenome: public ReadShiftProfile {
 private:
  int32_t numthreads;
  
 protected:
  int32_t ng_from;
  int32_t ng_to;
  int32_t ng_step;
  std::vector<range> seprange;
  
 public:
  std::string name;
  std::vector<ReadShiftProfile> chr;
  
 ReadShiftProfileGenome(std::string n, const Mapfile &p, const MyOpt::Variables &values):
  ReadShiftProfile(p.getlenF3(), values["ng_from"].as<int32_t>(), values["num4ssp"].as<int32_t>()),
    numthreads(values["threads"].as<int32_t>()),
    ng_from(5000),
    ng_to(values["ng_to"].as<int32_t>()),
    ng_step(values["ng_step"].as<int32_t>()),
    name(n)
    {
      for(auto x:p.genome.chr) {
	if(x.isautosome()) {
	  nread += x.bothnread_nonred();
	  len   += x.getlenmpbl();
	  //	  std::cout<< len << "\t" << x.getlenmpbl() << std::endl;
	}
      }
      for(auto x:p.genome.chr) {
	ReadShiftProfile v(p.getlenF3(), values["ng_from"].as<int32_t>(), values["num4ssp"].as<int32_t>(), 0, x.getlen(), x.bothnread_nonred(), x.getlenmpbl());
	v.setrchr(nread);
	chr.push_back(v);
      }
      // seprange
      defSepRange(numthreads);
    }
  virtual ~ReadShiftProfileGenome(){}
  void defSepRange(const int32_t numthreads) {
    int32_t length(mp_to+mp_from);
    int32_t sepsize = length/numthreads +1;
    for(int32_t i=0; i<numthreads; ++i) {
      int32_t s = i*sepsize;
      int32_t e = (i+1)*sepsize;
      if(i==numthreads-1) e = length;
      seprange.push_back(range(s - mp_from, e - mp_from));
    }
  }
  void addmp2genome(const int32_t i) {
    addmp(mp, chr[i].mp, chr[i].rchr);
    addmp(nc, chr[i].nc, chr[i].rchr);
  }

  void makeRscript(const std::string &filename, const std::string &prefix) {
    std::string Rscript(prefix + ".R");
    std::ofstream out(Rscript);

    out << "data <- read.csv('" << filename << "', header=TRUE, skip=6, sep='\t', quote='')" << std::endl;
    out << "output <- '" << prefix << "'" << std::endl;
    out << "pdf(paste(output, '.pdf', sep=''), height=7, width=14)" << std::endl;
    out << "par(mfrow=c(1,2))" << std::endl;
    out << "plot(data[1:" << mp_from+mp_to << ",1], data[1:" << mp_from+mp_to << ",5], type='l', xlab='Strand shift', ylab='Score relative to background', xlim=c(" << -mp_from << "," << mp_to << "), main='" << -mp_from << " bp ~ " << mp_to << " bp', ";
    if(name == "Jaccard index") out << "sub=sprintf('NSC=%g, RSC=%g, RLSC=%g, Bu=%g', " << nsc << "," << rsc << ","  << rlsc << "," << backgroundUniformity << "))" << std::endl;
    else if(name == "Cross correlation") out << "sub=sprintf('NSC=%g, RSC=%g, RLSC=%g', " << nsc << "," << rsc << ","  << rlsc << "))" << std::endl;
    else out << ")" << std::endl;
    out << "abline(v=" << nsci <<",lty=2,col=2)" << std::endl;
    out << "legend('bottomright', legend=paste('Estimated fragment length = ', " << nsci << "))" << std::endl;
    out << "plot(data[,1], data[,5], type='l', xlab='Strand shift',ylab='Score relative to background', main='Long distance (log scale)', log='x', xlim=c(1,1000000))" << std::endl;
    out << "dev.off()" << std::endl;

    std::string command = "R --vanilla < " + Rscript + " > " + Rscript + ".log 2>&1";
    
    int32_t return_code = system(command.c_str());
    if(WEXITSTATUS(return_code)) {
      std::cerr << "Warning: command " << command << "return nonzero status." << std::endl;
    }
  }
  
  void outputmpGenome(const std::string &prefix) {
    std::string filename = prefix + ".csv";
    print2file(filename, name);
    makeRscript(filename, prefix);
  }
  void outputmpChr(const std::string &filename, const int32_t i) {
    chr[i].print2file(filename, name);
  }
  void printStartMessage() const {
    std::cout << "Making " << name << " profile..." << std::flush;
  }
};

class shiftJacVec : public ReadShiftProfileGenome {
 public:
 shiftJacVec(const Mapfile &p, const MyOpt::Variables &values):
  ReadShiftProfileGenome("Jaccard index", p, values) {}

  void setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev);
  void execchr(const Mapfile &p, int32_t i) {
    auto fwd = genVector(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);  
  }
};

class shiftJacBit : public ReadShiftProfileGenome {
 public:
 shiftJacBit(const Mapfile &p, const MyOpt::Variables &values):
  ReadShiftProfileGenome("Jaccard index", p, values) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int32_t i) {
    auto fwd = genBitset(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftCcp : public ReadShiftProfileGenome {
 public:
 shiftCcp(const Mapfile &p, const MyOpt::Variables &values):
  ReadShiftProfileGenome("Cross correlation", p, values) {}

  void setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev);
  void execchr(const Mapfile &p, int32_t i) {
    auto fwd = genVector(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftHamming : public ReadShiftProfileGenome {
 public:
 shiftHamming(const Mapfile &p, const MyOpt::Variables &values):
  ReadShiftProfileGenome("Hamming distance", p, values) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int32_t i) {
    auto fwd = genBitset(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);
    
    setDist(chr[i], fwd, rev);
  }
};

class shiftFragVar : public ReadShiftProfileGenome {
  std::map<int32_t, PropNeighborFrag> distpnf;
  int32_t flen;
  double r4cmp;
  uint32_t numUsed4FCS;
  bool lackOfReads;
  int32_t ng_from_fcs;
  int32_t ng_to_fcs;
  int32_t ng_step_fcs;

 public:
 shiftFragVar(const Mapfile &p, const MyOpt::Variables &values, const int32_t fl):
  ReadShiftProfileGenome("Fragment Variability", p, values), flen(fl), r4cmp(0), numUsed4FCS(0), lackOfReads(false),
    ng_from_fcs(values["ng_from_fcs"].as<int32_t>()),
    ng_to_fcs(values["ng_to_fcs"].as<int32_t>()),
    ng_step_fcs(values["ng_step_fcs"].as<int32_t>())
      {
	//double r = (values["num4ssp"].as<int32_t>()/static_cast<double>(dist.getnread())) / (NUM_100M/static_cast<double>(dist.getlen()));
	double r = values["num4ssp"].as<int32_t>()/static_cast<double>(getnread());
#ifdef DEBUG
	std::cout << "\nr for FCS\t" << r << "\t reads: " << dist.getnread()<<  std::endl;
#endif
	if(r>1){
	  std::cerr << "\nWarning: number of reads (" << getnread() << ") is less than num4ssp ("<<  values["num4ssp"].as<int32_t>() <<").\n";
	  lackOfReads=true;
	}
	r4cmp = r*RAND_MAX;
      }
  std::vector<int8_t> genVector4FixedReadsNum(const strandData &seq, int32_t start, int32_t end);
  
  void setDist(ReadShiftProfile &chr, const std::vector<int8_t> &fwd, const std::vector<int8_t> &rev);
  void execchr(const Mapfile &p, const int32_t i) {
    auto fwd = genVector4FixedReadsNum(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector4FixedReadsNum(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
    addmp2genome(i);
  }

  uint32_t getnumUsed4FCS() const { return numUsed4FCS; }

  void printdistpnf(const std::string &filename) const {
    std::ofstream out(filename);

    for(auto x: v4pnf) out << "\tPNF len" << x;
    for(auto x: v4pnf) out << "\tCPNF len" << x;
    out << std::endl;

    for(size_t k=0; k<sizeOfvNeighborFrag-1; ++k) {
      out << k << "\t";
      for(auto x: v4pnf) out << distpnf.at(x).getPNF(k) << "\t";
      for(auto x: v4pnf) out << distpnf.at(x).getCumulativePNF(k) << "\t";
      out << std::endl;
    }
  }
  
  void outputfcsGenome(const std::string &filename) const {
    print2file4fcs(filename, name, flen, lackOfReads);
  }
  void outputfcsChr(const std::string &filename, const int32_t i) const {  
    chr[i].print2file4fcs(filename, name, flen, lackOfReads);
  }

  double getMPflen() const { return mp.at(flen); } 
};

#endif /* _SSP_SHIFTPROFILE_P_H_ */
