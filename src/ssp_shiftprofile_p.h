/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#ifndef _SSP_SHIFTPROFILE_P_H_
#define _SSP_SHIFTPROFILE_P_H_

#include "pw_gv.h"
#include <boost/dynamic_bitset.hpp>
#include "alglib/alglib.h"

namespace {
  const int mp_from(500);
  const int mp_to(1500);
  const int sizeOfvDistOfDistaneOfFrag = 5000;
  const std::vector<int> v4acfp{50, 100, 150, 500, 1000, 2000, 3000, 10000, 100000, 1000000};
}

std::vector<char> genVector(const strandData &seq, int start, int end);
std::vector<char> genVector4FixedReadsNum(const strandData &seq, int start, int end, const double r4cmp);
boost::dynamic_bitset<> genBitset(const strandData &seq, int, int);
void addmp(std::map<int, double> &, const std::map<int, double> &, double w);

class FragmentVariability {
  double sumOfvDistOfDistaneOfFrag;
  std::vector<int> vDistOfDistaneOfFrag;

 public:
 FragmentVariability(): sumOfvDistOfDistaneOfFrag(0),
    vDistOfDistaneOfFrag(sizeOfvDistOfDistaneOfFrag, 0) {}

  void setVariability(const int fraglen, const int start, const int end,
		      const std::vector<char> &fwd, const std::vector<char> &rev) {
    int s(start);
    if(start + fraglen < 0) s = start - fraglen;
    int e(end - fraglen);
    int last(s);
    for(int i=s; i<e; ++i) {
      if(fwd[i] && rev[i+fraglen]) {
	if(RANGE(i-last, 0, sizeOfvDistOfDistaneOfFrag-1)) ++vDistOfDistaneOfFrag[i-last];
	else ++vDistOfDistaneOfFrag[sizeOfvDistOfDistaneOfFrag-1];
	last = i;
      }
    }

    sumOfvDistOfDistaneOfFrag = accumulate(vDistOfDistaneOfFrag.begin(), vDistOfDistaneOfFrag.end(), 0);
  }

  double getsumOfvDistOfDistaneOfFrag() const { return sumOfvDistOfDistaneOfFrag; }
  double getDistOfDistanceOfFragment(const int i) const {
    if(i<0 || i>= sizeOfvDistOfDistaneOfFrag) {
      std::cerr << "error: invalid num " << i << "for getDistOfDistanceOfFragment max: " << sizeOfvDistOfDistaneOfFrag << std::endl;
      return -1;
    }
    else {
      return sumOfvDistOfDistaneOfFrag ? vDistOfDistaneOfFrag[i] / sumOfvDistOfDistaneOfFrag : 0;
    }
  }
  double getAccuOfDistanceOfFragment(const int i) const {
    if(i<0 || i>= sizeOfvDistOfDistaneOfFrag) {
      std::cerr << "error: invalid num " << i << "for getAccuOfDistanceOfFragment max: " << sizeOfvDistOfDistaneOfFrag << std::endl;
      return -1;
    }
    else {
      double vAccuOfDistaneOfFrag(0);
      for(int j=0; j<=i; ++j) {
	vAccuOfDistaneOfFrag += getDistOfDistanceOfFragment(j);
      }
      return vAccuOfDistaneOfFrag;
    }
  }
  double getPvalueOfBinomTest(const int i, const double myu) const {
    if(i<0 || i>= sizeOfvDistOfDistaneOfFrag) {
      std::cerr << "error: invalid num " << i << "for getPvalueOfBinomTest max: " << sizeOfvDistOfDistaneOfFrag << std::endl;
      return -1;
    }
    else {
      double sum(0);
      for(int j=0; j<=i; ++j) {
	sum += vDistOfDistaneOfFrag[j];
      }
      std::cout << sum << "," <<  sumOfvDistOfDistaneOfFrag << "," << myu << std::endl;
      return stdNormdist(sum, sumOfvDistOfDistaneOfFrag*myu, sqrt(sumOfvDistOfDistaneOfFrag*myu*(1-myu))); 
    }
  }
  void add2genome(const FragmentVariability &x, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    for(uint i=0; i<vDistOfDistaneOfFrag.size(); ++i) vDistOfDistaneOfFrag[i] += x.vDistOfDistaneOfFrag[i];
    sumOfvDistOfDistaneOfFrag += x.sumOfvDistOfDistaneOfFrag;
  }
};

class ReadShiftProfile {
  int lenF3;
  double r;
  double bk;
  int bk_from;

 protected:
  double nsc;
  double rlsc;
  int nsci;
  long len;
  long nread;
  double backgroundUniformity;
  
 public:
  std::map<int, double> mp;
  std::map<int, double> nc;
  int start;
  int end;
  int width;

  double rchr;

 ReadShiftProfile(const int len, const int b, int s=0, int e=0, long n=0, long l=0):
  lenF3(len), r(0), bk(0), bk_from(b), nsc(0), rlsc(0), nsci(0), len(l), nread(n), backgroundUniformity(0), start(s), end(e), width(e-s), rchr(1) {}
  virtual ~ReadShiftProfile() {}
  void setmp(const int i, const double val, boost::mutex &mtx) {
    boost::mutex::scoped_lock lock(mtx);
    mp[i] = val;
  }

  double getbackgroundUniformity() const { return backgroundUniformity; }
  void setrchr(const long n) { rchr = n ? nread/static_cast<double>(n): 0; }
  int getlenF3() const { return lenF3; }
  int getnsci() const { return nsci; }
  double getmpsum() const {
    double sum(0);
    for(auto itr = mp.begin(); itr != mp.end(); ++itr) sum += itr->second;
    return sum;
  }
  void setControlRatio() {
    int n(0);
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
    int threwidth(5);
    
    setControlRatio();
    nsc = mp[mp_to-1];
    for(int i=mp_to-1-threwidth; i > lenF3*1.3; --i) {
      int on(1);
      if(name == "Hamming distance") {
	for(int j=1; j<=threwidth; ++j) {
	  if (mp[i] > mp[i+j] || mp[i] > mp[i-j]) on=0;
	}
	if(on && nsc > mp[i]) {
	  nsc  = mp[i];
	  nsci = i;
	}
      } else {
	for(int j=1; j<=threwidth; ++j) {
	  if (mp[i] < mp[i+j] || mp[i] < mp[i-j]) on=0;
	}
	if(on && nsc < mp[i] * r) {
	  nsc  = mp[i] * r;
	  nsci = i;
	}
      }
    }
  }

  void print2file(const std::string &filename, const std::string &name) {
    if(!nread) {
      std::cerr << filename << ": no read" << std::endl;
    }
    double sum(getmpsum());


    ////////////////
    double rRPKM = (NUM_10M/static_cast<double>(nread)) / (NUM_100M/static_cast<double>(len));
    double be(bk * rRPKM);
    double const_bu(1/39.0);  // N/(4*L-N), N=10M, L=100M
    rlsc = mp.at(lenF3) *r;
    backgroundUniformity = const_bu / be;
    
    std::ofstream out(filename);
    out << "NSC\t" << nsc  << std::endl;
    out << "RLSC\t"<< rlsc << std::endl;
    out << "Estimated fragment length\t" << nsci << std::endl;
    out << "Background enrichment\t" << be << std::endl;
    out << "Background uniformity\t" << backgroundUniformity << std::endl;

    out << "Strand shift\t" << name << "\tprop\tper 10M reads\tper control" << std::endl;
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
  void print2file4fcs(const std::string filename, const std::string name, const int flen, const bool lackOfReads) const {
    if(!nread) {
      std::cerr << filename << ": no read" << std::endl;
    }
    std::ofstream out(filename);

    //    std::cout << flen <<"," << lenF3 <<std::endl;

    std::string str("");
    if(lackOfReads) str = " (read number is insufficient)";
    out << "Fragment score" << str << "\t" << mp.at(flen) << std::endl;
    out << "Read score" << str << "\t" << mp.at(lenF3) << std::endl;
    out << "Strand shift\t" << name << std::endl;
    for(auto itr = mp.begin(); itr != mp.end(); ++itr)
      if(itr->first != flen && itr->first != lenF3)
	out << itr->first << "\t" << itr->second << std::endl;
  }
};

class ReadShiftProfileGenome: public ReadShiftProfile {
 private:
  int numthreads;
  
 protected:
  int ng_from;
  int ng_to;
  int ng_step;
  std::vector<range> seprange;
  
 public:
  std::string name;
  std::vector<ReadShiftProfile> chr;
  
 ReadShiftProfileGenome(std::string n, const Mapfile &p, const MyOpt::Variables &values):
  ReadShiftProfile(p.getlenF3(), values["ng_from"].as<int>()),
    numthreads(values["threads"].as<int>()),
    ng_from(5000),
    ng_to(values["ng_to"].as<int>()),
    ng_step(values["ng_step"].as<int>()),
    name(n) {
    for(auto x:p.genome.chr) {
      if(x.isautosome()) {
	nread += x.bothnread_nonred();
	len   += x.getlenmpbl();
      }
    }
    for(auto x:p.genome.chr) {
      ReadShiftProfile v(p.getlenF3(), values["ng_from"].as<int>(), 0, x.getlen(), x.bothnread_nonred(), x.getlenmpbl());
      v.setrchr(nread);
      chr.push_back(v);
    }
    // seprange
    defSepRange(numthreads);
  }
  virtual ~ReadShiftProfileGenome(){}
  long getnread() const { return nread; }
  void defSepRange(const int numthreads) {
    int length(mp_to+mp_from);
    int sepsize = length/numthreads +1;
    for(int i=0; i<numthreads; ++i) {
      int s = i*sepsize;
      int e = (i+1)*sepsize;
      if(i==numthreads-1) e = length;
      seprange.push_back(range(s - mp_from, e - mp_from));
    }
  }
  void addmp2genome(const int i) {
    addmp(mp, chr[i].mp, chr[i].rchr);
    addmp(nc, chr[i].nc, chr[i].rchr);
  }

  void makeRscript(const std::string &filename, const std::string &prefix) {
    std::string Rscript(prefix + ".R");
    std::ofstream out(Rscript);

    out << "data <- read.csv('" << filename << "', header=TRUE, skip=5, sep='\t', quote='')" << std::endl;
    out << "output <- '" << prefix << "'" << std::endl;
    out << "pdf(paste(output, '.pdf', sep=''), height=7, width=14)" << std::endl;
    out << "par(mfrow=c(1,2))" << std::endl;
    if(name == "Jaccard index") out << "plot(data[,1], data[,4], type='l', xlab='Strand shift', ylab='Normalized score', xlim=c(" << -mp_from << "," << mp_to << "), main='" << -mp_from << " bp ~ " << mp_to << " bp', sub=sprintf('NSC=%g, RLSC=%g, Bu=%g', " << nsc << "," << rlsc << "," << backgroundUniformity << "))" << std::endl;
    else if(name == "Cross correlation") out << "plot(data[,1], data[,4], type='l', xlab='Strand shift', ylab='Normalized score', xlim=c(-200,1500), main='-200 bp ~ 1500 bp', sub=sprintf('NSC=%g, RLSC=%g', " << nsc << "," << rlsc << "))" << std::endl;
    else out << "plot(data[,1], data[,4], type='l', xlab='Strand shift', ylab='Normalized score', xlim=c(-200,1500), main='-200 bp ~ 1500 bp')" << std::endl;
    out << "abline(v=" << nsci <<",lty=2,col=2)" << std::endl;
    out << "legend('bottomright', legend=paste('Estimated fragment length = ', " << nsci << "))" << std::endl;
    out << "plot(data[,1], data[,4], type='l', xlab='Strand shift',ylab='Normalized score', main='Long distance')" << std::endl;
    out << "dev.off()" << std::endl;

    std::string command = "R --vanilla < " + Rscript + " | tee " + Rscript + ".log";
    
    int return_code = system(command.c_str());
  }
  
  void outputmpGenome(const std::string &prefix) {
    std::string filename = prefix + ".csv";
    print2file(filename, name);
    makeRscript(filename, prefix);
  }
  void outputmpChr(const std::string &filename, const int i) {
    chr[i].print2file(filename, name);
  }
  void printStartMessage() const {
    std::cout << "Making " << name << " profile..." << std::flush;
  }
  
  void printacfp(const std::string &){};
};

class shiftJacVec : public ReadShiftProfileGenome {
 public:
 shiftJacVec(const Mapfile &p, const MyOpt::Variables &values): ReadShiftProfileGenome("Jaccard index", p, values) {}

  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);  
  }
};

class shiftJacBit : public ReadShiftProfileGenome {
 public:
 shiftJacBit(const Mapfile &p, const MyOpt::Variables &values): ReadShiftProfileGenome("Jaccard index", p, values) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftCcp : public ReadShiftProfileGenome {
 public:
 shiftCcp(const Mapfile &p, const MyOpt::Variables &values): ReadShiftProfileGenome("Cross correlation", p, values) {}

  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genVector(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genVector(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);

    setDist(chr[i], fwd, rev);
  }
};

class shiftHamming : public ReadShiftProfileGenome {
 public:
 shiftHamming(const Mapfile &p, const MyOpt::Variables &values): ReadShiftProfileGenome("Hamming distance", p, values) {}

  void setDist(ReadShiftProfile &chr, const boost::dynamic_bitset<> &fwd, boost::dynamic_bitset<> &rev);
  void execchr(const Mapfile &p, int i) {
    auto fwd = genBitset(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end);
    auto rev = genBitset(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end);
    
    setDist(chr[i], fwd, rev);
  }
};

class shiftFragVar : public ReadShiftProfileGenome {
  std::map<int, FragmentVariability> acfp;
  int flen;
  bool lackOfReads;
  bool fcsfull;
  int ng_from_fcs;
  int ng_to_fcs;
  int ng_step_fcs;

 public:
 shiftFragVar(const Mapfile &p, const MyOpt::Variables &values, const int fl):
  ReadShiftProfileGenome("Fragment Variability", p, values), flen(fl), lackOfReads(false), fcsfull(values.count("fcsfull")),
    ng_from_fcs(values["ng_from_fcs"].as<int>()),
    ng_to_fcs(values["ng_to_fcs"].as<int>()),
    ng_step_fcs(values["ng_step_fcs"].as<int>()) {}

  void setDist(ReadShiftProfile &chr, const std::vector<char> &fwd, const std::vector<char> &rev);
  void execchr(const Mapfile &p, const int i, const double r4cmp) {
    auto fwd = genVector4FixedReadsNum(p.genome.chr[i].seq[STRAND_PLUS],  chr[i].start, chr[i].end, r4cmp);
    auto rev = genVector4FixedReadsNum(p.genome.chr[i].seq[STRAND_MINUS], chr[i].start, chr[i].end, r4cmp);

    setDist(chr[i], fwd, rev);
  }

  void lackOfReads_on() { lackOfReads=true; }
  void printacfp(const std::string &filename) const {
    std::ofstream out(filename);

    for(auto x: v4acfp) {
      if(fcsfull && x > mp_to) continue;
      out << "\tlen" << x;
    }
    out << std::endl;

    /*    for(auto x: v4acfp) {
      double myu(acfp.at(100000).getAccuOfDistanceOfFragment(sizeOfvDistOfDistaneOfFrag-2));
      if(fcsfull && x > mp_to) continue;
      out << "\t" << acfp.at(x).getPvalueOfBinomTest(sizeOfvDistOfDistaneOfFrag-2, myu);
    }
    out << std::endl;*/

    for(size_t k=0; k<sizeOfvDistOfDistaneOfFrag-1; ++k) {
      out << k << "\t";

      for(auto x: v4acfp) {
	if(fcsfull && x > mp_to) continue;
	out << acfp.at(x).getAccuOfDistanceOfFragment(k) << "\t";
      }
      out << std::endl;
    }
  }
  
  void outputmpGenome(const std::string &filename) const {
    print2file4fcs(filename, name, flen, lackOfReads);
  }
  void outputmpChr(const std::string &filename, const int i) const {  
    chr[i].print2file4fcs(filename, name, flen, lackOfReads);
  }
};

#endif /* _SSP_SHIFTPROFILE_P_H_ */
