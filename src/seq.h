/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef SEQ_H
#define SEQ_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

std::string rmchr(const std::string &chr);

enum Strand {STRAND_PLUS, STRAND_MINUS, STRANDNUM, STRAND_BOTH};
enum status {INTERGENIC, GENIC, INTRON, EXON, DOWNSTREAM, UPSTREAM, TSS, PARALLEL, DIVERGENT, CONVERGENT};

template <class T>
class var {
  std::string name;
  T val;
  T limlow;
  T limup;
  bool isupper;
 public:
 var(): name(""), val(0), limlow(0), limup(0){}
  var(std::string str, T low):       name(str), val(0), limlow(low), limup(0), isupper(false) {}
  var(std::string str, T low, T up): name(str), val(0), limlow(low), limup(up), isupper(true) {}
  void set(T n) {
    if(isupper && (n<limlow || n>limup)) {
      std::cout << "Error : variable " << name << " should be " << limlow << "<= and <=" << limup << "." << std::endl;
    }else if(!isupper && n<limlow ) {
      std::cout << "Error : variable " << name << " should be >=" << limlow << "." << std::endl;
    }
    else val=n;
  }
  operator T() const { return val; }
};

class chrsize {
  std::string name;
  int32_t len;

 public:
 chrsize(): name(""), len(0) {}
 chrsize(const std::string &n, const int32_t l): name(n), len(l) {}
  std::string getname() const { return name; }
  int32_t getlen() const { return len; }
};

class sepchr {
 public:
  uint32_t s;
  uint32_t e;
 sepchr(uint32_t start, uint32_t end): s(start), e(end) {}
};

class range {
 public:
  int32_t start;
  int32_t end;
 range(): start(0), end(0) {}
 range(int32_t s, int32_t e): start(s), end(e) {}
};

class bed {
 public:
  std::string chr;
  int32_t start;
  int32_t end;
  int32_t summit;
 bed(): start(0), end(0), summit(0) {}
  virtual ~bed(){}
 bed(int32_t s, int32_t e, std::string c): chr(c), start(s), end(e) {}
 bed(std::vector<std::string> s): start(stoi(s[1])), end(stoi(s[2])), summit((start + end)/2) {
    if(!s[0].find("chr")) chr = s[0].substr(3);
    else chr = s[0];
  }
  void print()      const { std::cout << "chr" << chr << "\t" << start  << "\t" << end ; }
  void printHead () const { std::cout << "chromosome\tstart\tend"; }
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
 void printHead () const {
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

class genedata {
 public:
  std::string tname;
  std::string gname;
  std::string tid;
  std::string gid;
  std::string chr;
  int32_t txStart;   // "Transcription start position"
  int32_t txEnd;     // "Transcription end position"
  int32_t cdsStart;  // "Coding region start"
  int32_t cdsEnd;    // "Coding region end"
  int32_t exonCount; // "Number of exons"
  std::string strand;
  std::vector<range> exon;

  // for Ensembl
  std::string gsrc;  // gene source
  std::string tsrc;  // transcript source
  std::string gtype; // gene biotype
  std::string ttype; // transcript biotype
  std::string ttag;  // Gencode tag

  genedata(): txStart(0), txEnd(0), cdsStart(0), cdsEnd(0), exonCount(0) {}

  int32_t length() const { return (txEnd - txStart); }
  void printall() const {
    if(this){
      std::cout << tname << "\t" << gname << "\t" << tid << "\t" << gid << "\t" << chr << "\t" << strand << "\t" << txStart << "\t" << txEnd << "\t" << cdsStart << "\t" << cdsEnd << "\t" << exonCount << "\tgene source: " << gsrc << "\ttranscript source: "<< tsrc << "\tgene biotype: "<< gtype << "\ttranscript biotype: "<< ttype  << "\ttranscript tag: "<< ttag << "\t";
      for (auto x: exon) std::cout << x.start << "-" << x.end << ", ";
    }
  }
  void print() const {
    if(this) std::cout << tname << "\t" << gname << "\t" << tid << "\t" << gid << "\t" << strand << "\t" << txStart << "\t" << txEnd << "\t";
  }
};

class hitGene {
 public:
  status st;
  int32_t d;
  const genedata *gene;
 hitGene(): st(INTERGENIC), d(0), gene(nullptr) {}
};

template <class T>
class bed_gene {
 public:
  T bed;
  hitGene gene;
  std::vector<hitGene> genelist;
 bed_gene(): gene() {}
 bed_gene(std::vector<std::string> s): bed(s), gene() {}
  void print() const { bed.print();}
  void printWithGene(bool redundant) const {
    if(redundant) {
      for(auto x:genelist) {
	print();
	if(x.st == UPSTREAM)        std::cout << "\tupstream\t";
	else if(x.st == DOWNSTREAM) std::cout << "\tdownstream\t";
	else if(x.st == GENIC)      std::cout << "\tgenic\t";
	else if(x.st == INTERGENIC) std::cout << "\tintergenic\t";
	else if(x.st == CONVERGENT) std::cout << "\tconvergent\t";
	else if(x.st == DIVERGENT)  std::cout << "\tdivergent\t";
	else if(x.st == PARALLEL)   std::cout << "\tparallel\t";
	x.gene->print();
	std::cout << std::endl;
      }
    } else {
      print();
      if(gene.st == UPSTREAM)        std::cout << "\tupstream\t";
      else if(gene.st == DOWNSTREAM) std::cout << "\tdownstream\t";
      else if(gene.st == GENIC)      std::cout << "\tgenic\t";
      else if(gene.st == INTERGENIC) std::cout << "\tintergenic\t";
      else if(gene.st == CONVERGENT) std::cout << "\tconvergent\t";
      else if(gene.st == DIVERGENT)  std::cout << "\tdivergent\t";
      else if(gene.st == PARALLEL)   std::cout << "\tparallel\t";
      gene.gene->print();
      std::cout << std::endl;
    }
  }
  void printWithTss(bool redundant) const {
    if(redundant) {
      for(auto x:genelist) {
	print();
	std::cout << "\t" << x.d << "\t"; 	
	x.gene->print();
	std::cout << std::endl;
      }
    } else {
      print();
      if(gene.st == TSS) {
	std::cout << "\t" << gene.d << "\t"; 	
	gene.gene->print();
      }
      std::cout << std::endl;
    }
  }
  
  void update(const status &pst, const genedata &pgene) {
    if(gene.st < pst){
      gene.st = pst;
      gene.gene = &pgene;
    }
    hitGene g;
    g.st = pst;
    g.gene = &pgene;
    genelist.push_back(g);
  }
  void update(const status &pst, const genedata *pgene) {
    if(gene.st < pst){
      gene.st = pst;
      gene.gene = pgene;
    }
    hitGene g;
    g.st = pst;
    g.gene = pgene;
    genelist.push_back(g);
  }
  void printHead () const { bed.printHead(); }
};

class strRange : public bed {
 public:
  std::string strand;
  const genedata *gene;
 strRange(): bed(), strand(0), gene(0) {}
 strRange(int32_t s, int32_t e, std::string c, std::string str, const genedata &pgene): bed(s,e,c), strand(str), gene(&pgene) {}
};

class convsite : public bed {
 public:
  status st;
  const genedata *gene;
 convsite(int32_t s, int32_t e, std::string c, const genedata *pgene): bed(s,e,c), gene(pgene) {}
};

class fasta {
 public:
  std::string name;
  uint64_t len, len_mpbl;
  int32_t nbin;
  double p_mpbl;  /* mappability */
  double gcov;    /* genome coverage for bin */
  fasta (std::string str, int32_t l=0): name(str), len(l), len_mpbl(0), nbin(0), p_mpbl(0), gcov(0) {}
  fasta (std::vector<std::string> &v): name(v[0]), len(stoi(v[1])), len_mpbl(0), p_mpbl(0), gcov(0) {}
  void print() const {
    std::cout << name << "\t" << len << "\t" << nbin << "\t" << len_mpbl << "\t"<< p_mpbl << "\t" << gcov << std::endl;
  }
};

class Peak : public bed {
 public:
  int32_t summit;
  double pileup;
  double enrich;
  double p_inter, p_enr;
  double q;
 Peak(int32_t s, int32_t e, std::string c, double val, double p):
  bed(s,e,c), summit(s), pileup(val), enrich(0), p_inter(p), p_enr(0), q(0) {}
  void renew(int32_t i, double val, double p) {
    end = i;
    pileup += val;
    if(p_inter > p) {
      p_inter = p;
      summit = i;
    }
  }
  void print(std::ofstream &out, int32_t id, int32_t binsize) const {
    out << chr << "\t" << start*binsize << "\t" << end*binsize << "\t"
	<< ((end - start +1)*binsize-1) << "\t" << (summit*binsize -binsize/2) << "\t" << pileup << "\t"
	<< p_inter << "\t" << enrich << "\t" << q << "\tpeak " << id << std::endl;
  }
  void printHead (std::ofstream &out) const {
    out << "chromosome\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname" << std::endl;
  }
};


class Fragment {
public:
  std::string chr;
  int32_t F3;
  Strand strand;
  int32_t fraglen;
  int32_t readlen_F3;

 Fragment(): fraglen(0), readlen_F3(0) {}
 void addSAM(const std::vector<std::string> &v, const bool pair, const int32_t sv) {
   chr = rmchr(v[2]);
   readlen_F3 = v[9].length();
   if(pair) fraglen = abs(stoi(v[8]));
   else fraglen = readlen_F3;
   if(sv&16) {
     strand = STRAND_MINUS;
     F3 = stoi(v[3]) + readlen_F3 -1;
   } else {
     strand = STRAND_PLUS;
     F3 = stoi(v[3]) -1;
   }
 }
 void print() const {
   std::cout << "chr:" << chr
	     << "\tposi:" << F3
	     << "\tstrand:" << strand
	     << "\tfraglen:" << fraglen
	     <<"\treadlen:" << readlen_F3
	     << std::endl;
  }
};

class Read {
  int32_t weight;
  enum {WeightNum=1000};
 public:
  int32_t F3;
  int32_t F5;
  int32_t duplicate;
  int32_t inpeak;
  
 Read(const Fragment &frag): weight(WeightNum), F3(frag.F3), duplicate(0), inpeak(0) {
    if(frag.strand == STRAND_PLUS) F5 = frag.F3 + frag.fraglen;
    else F5 = frag.F3 - frag.fraglen;
  }
  double getWeight() const {
    return weight/static_cast<double>(WeightNum);
  }
  void multiplyWeight(const double w) { weight *= w; }
};


#endif  // SEQ_H
