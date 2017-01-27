/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <ext/stdio_filebuf.h>
#include "pw_readmapfile.h"

namespace {
  void addFragToChr(SeqStatsGenome &genome, const Fragment &frag)
  {
    genome.dflen.addF3(frag.readlen_F3);
    genome.dflen.addvflen(frag.fraglen);
    
    bool on(false);
    for(auto &x: genome.chr) {
      if(x.getname() == frag.chr) {
	x.addfrag(frag);
	on=true;
      }
    }
    if(!on) std::cerr << "Warning: " << frag.chr << " is not in genometable." << std::endl;
    
    return;
  } 
  
  template <class T>
  void do_bampe(const MyOpt::Variables &values, SeqStatsGenome &genome, T &in)
  {
    int32_t maxins(values["maxins"].as<int32_t>());
    
    std::string lineStr;
    while (!in.eof()) {
      getline(in, lineStr);
      if(lineStr.empty() || lineStr[0]=='@') continue;
      
      std::vector<std::string> v;
      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      int32_t sv(stoi(v[1]));   // bitwise FLAG
      if(sv&4 || sv&512 || sv&1024) continue;
      if(!(sv&2)) continue;
      if(sv&128) {
	genome.dflen.addF5(v[9].length());
	continue;
      }
      Fragment frag;
      frag.addSAM(v, values.count("pair"), sv);
      if(frag.fraglen > maxins) continue;
      //frag.print();
      
      addFragToChr(genome, frag);
    }
    
    return;
  }
  
  template <class T>
  void do_bamse(const MyOpt::Variables &values, SeqStatsGenome &genome, T & in)
  {
    std::string lineStr; 
    while (!in.eof()) {
      getline(in, lineStr);
      if(lineStr.empty() || lineStr[0]=='@') continue;
      std::vector<std::string> v;
      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      int32_t sv(stoi(v[1])); // bitwise FLAG
      // unmapped reads, low quality reads
      if(sv&4 || sv&512 || sv&1024) continue;
      if(sv&64 || sv&128) std::cerr << "Warning: parsing paired-end file as single-end." << std::endl;
      Fragment frag;
      frag.addSAM(v, values.count("pair"), sv);
      //    std::cout << lineStr << std::endl;
      // frag.print();
      addFragToChr(genome, frag);
    }
    
    return;
  }
  
  void parseSam(const MyOpt::Variables &values, const std::string &inputfile, SeqStatsGenome &genome)
  {
    if((values.count("ftype") && values["ftype"].as<std::string>()=="SAM") || isStr(inputfile, ".sam")) {  // SAM
      std::cout << "Input format: SAM" << std::endl;
      std::ifstream in(inputfile);
      if(!in) PRINTERR("Could not open " << inputfile << ".");
      if (values.count("pair")) do_bampe(values, genome, in);
      else do_bamse(values, genome, in);
    }
    else {  // BAM
      std::cout << "Input format: BAM" << std::endl;
      std::string command = "samtools view -h " + inputfile;
      FILE *fp = popen(command.c_str(), "r");
      __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
      std::istream in(static_cast<std::streambuf *>(p_fb));
      if (values.count("pair")) do_bampe(values, genome, in);
      else do_bamse(values, genome, in);
    }
    
    return;
  }
  
  void parseBowtie(const MyOpt::Variables &values, const std::string &inputfile, SeqStatsGenome &genome)
  {
    int32_t maxins(values["maxins"].as<int32_t>());
    std::ifstream in(inputfile);
    if(!in) PRINTERR("Could not open " << inputfile << ".");
    std::cout << "Input format: BOWTIE" << std::endl;
    
    std::string chr_F3(""), chr_F5(""), nametemp("");
    int32_t F5(0);
    Fragment fragpair;
    
    std::string lineStr;
    while (!in.eof()) {
      getline(in, lineStr);
      if(lineStr.empty()) continue;
      
      std::vector<std::string> v;
      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      
      if (values.count("pair")) {
	std::vector<std::string> read;
	boost::split(read, v[0], boost::algorithm::is_any_of("/"));
	if(nametemp != "" && nametemp != read[0]) PRINTERR("Invalid read pair." << nametemp <<"-" << read[0]);
	if(read[1] == "1") {  // F3 read
	  chr_F3 = rmchr(v[2]);
	  fragpair.readlen_F3 = v[4].length();
	  if(v[1] == "+") { 
	    fragpair.strand = Strand::FWD;
	    fragpair.F3 = stoi(v[3]);
	  } else {
	    fragpair.strand = Strand::REV;
	    fragpair.F3 = stoi(v[3]) + fragpair.readlen_F3;
	  }
	} else {  
	  chr_F5 = rmchr(v[2]);
	  if(v[1] == "+") F5 = stoi(v[3]);
	  else            F5 = stoi(v[3]) + v[4].length();
	  genome.dflen.addF5(v[4].length());
	}
	if(chr_F3 != "" && chr_F5 != ""){
	  if(chr_F3 == chr_F5) {
	    fragpair.chr = chr_F3;
	    fragpair.fraglen = abs(F5 - fragpair.F3);
	    if(fragpair.fraglen <= maxins) addFragToChr(genome, fragpair); 
	    //	  fragpair.print();
	  }
	  chr_F3 = "";
	  chr_F5 = "";
	  nametemp = "";
	}
      } else {
	if(isStr(v[0], "/2")) PRINTERR("Warning: parsing paired-end file as single-end");
	Fragment frag;
	frag.chr = rmchr(v[2]);
	frag.readlen_F3 = v[4].length();
	if(v[1] == "+") { 
	  frag.strand = Strand::FWD;
	  frag.F3 = stoi(v[3]);
	} else {
	  frag.strand = Strand::REV;
	  frag.F3 = stoi(v[3]) + frag.readlen_F3;
	}
	//      frag.print();
	addFragToChr(genome, frag); 
      }
      
    }
    return;
  }
  
  template <class T>
  void funcTagAlign(const MyOpt::Variables &values, SeqStatsGenome &genome, T &in)
  {
    std::cout << "Input format: TAGALIGN" << std::endl;
    std::string lineStr;
    while (!in.eof()) {
      getline(in, lineStr);
      if(lineStr.empty()) continue;
      
      std::vector<std::string> v;
      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      if(v.size() < 6) PRINTERR("Use tagAlign (BED3+3) file");
      
      if (values.count("pair")) PRINTERR("tagAlign format does not support paired-end file.\n");
      else {
	int32_t start(stoi(v[1]));
	int32_t end(stoi(v[2]));
	Fragment frag;
	frag.chr = rmchr(v[0]);
	frag.readlen_F3 = abs(end - start);
	if(v[5] == "+") {
	  frag.strand = Strand::FWD;
	  frag.F3 = start;
	} else {
	  frag.strand = Strand::REV;
	  frag.F3 = start + frag.readlen_F3;
	}
	//      std::cout << lineStr << std::endl;
	//frag.print();
	addFragToChr(genome, frag);
      }
    }
    return;
  }
  
  void parseTagAlign(const MyOpt::Variables &values, const std::string &inputfile, SeqStatsGenome &genome)
  {
    if(isStr(inputfile, ".gz")) {
      std::string command = "zcat " + inputfile;
      FILE *fp = popen(command.c_str(), "r");
      __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
      std::istream in(static_cast<std::streambuf *>(p_fb));
      funcTagAlign(values, genome, in);
    } else {
      std::ifstream in(inputfile);
      if(!in) PRINTERR("Could not open " << inputfile << ".");
      funcTagAlign(values, genome, in);
    }
    return;
  }

  void hashFilterAllSingle(std::unordered_map<int32_t, int32_t> &mp, strandData &seq, const int32_t thre)
  {
    for(auto &x:seq.vRead) {
      if(mp.find(x.F3) != mp.end()) {
	if(mp[x.F3] < thre) {
	  ++mp[x.F3];
	  seq.nread_nonred++;
	} else {
	  x.duplicate = 1;
	  seq.nread_red++;
	}
      } else {
	mp[x.F3] = 1;
	seq.nread_nonred++;
      }
    }
    return;
  }
  
  void hashFilterCmpSingle(std::unordered_map<int32_t, int32_t> &mp, class LibComp &comp, const strandData &seq, const int32_t thre)
  {
    for(auto &x: seq.vRead){
      if(rand() >= comp.getr4cmp()) continue;
      comp.incNtAll();
      if(mp.find(x.F3) != mp.end()) {
	if(mp[x.F3] < thre) {
	  ++mp[x.F3];
	  comp.incNtNonred();
	} else {
	  comp.incNtRed();
	}
      } else {
	mp[x.F3] = 1;
	comp.incNtNonred();
      }
    }
    return;
  }
  
  void hashFilterAllPair(std::unordered_map<std::string, int32_t> &mp, strandData &seq, const int32_t thre)
  {
    for(auto &x:seq.vRead) {
      int32_t Fmin = std::min(x.F3, x.F5);
      int32_t Fmax = std::max(x.F3, x.F5);
      std::string str = IntToString(Fmin) + "-" + IntToString(Fmax);
      //    std::cout << str << std::endl;
      if(mp.find(str) != mp.end()) {
	if(mp[str] < thre) {
	  ++mp[str];
	  ++seq.nread_nonred;
	} else {
	  x.duplicate = 1;
	  ++seq.nread_red;
	}
      } else {
	mp[str] = 1;
	++seq.nread_nonred;
      }
    }
    return;
  }
  
  void hashFilterCmpPair(std::unordered_map<std::string, int32_t> &mp, class LibComp &comp, const strandData &seq, const int32_t thre)
  {
    for(auto &x: seq.vRead){
      if(rand() >= comp.getr4cmp()) continue;
      int32_t Fmin = std::min(x.F3, x.F5);
      int32_t Fmax = std::max(x.F3, x.F5);
      std::string str = IntToString(Fmin) + "-" + IntToString(Fmax);
      comp.incNtAll();
      if(mp.find(str) != mp.end()) {
	if(mp[str] < thre) {
	  ++mp[str];
	  comp.incNtNonred();
	} else {
	  comp.incNtRed();
	}
      } else {
	mp[str] = 1;
	comp.incNtNonred();
      }
    }
    return;
  }
  
  void filtering_eachchr_single(Mapfile &p, SeqStats &chr)
  {
    for (auto strand: {Strand::FWD, Strand::REV}) {
      std::unordered_map<int32_t, int32_t> mp;
      hashFilterAllSingle(mp, chr.getStrandref(strand), p.complexity.getThreshold());
      
      std::unordered_map<int32_t, int32_t> mp2;
      hashFilterCmpSingle(mp2, p.complexity, chr.getStrandref(strand), p.complexity.getThreshold());
    }
    
    return;
  }
  
  void filtering_eachchr_pair(Mapfile &p, SeqStats &chr)
  {
    std::unordered_map<std::string, int32_t> mp;
    
    for (auto strand: {Strand::FWD, Strand::REV}) {
      hashFilterAllPair(mp, chr.getStrandref(strand), p.complexity.getThreshold());
    }
    
    std::unordered_map<std::string, int> mp2;
    for (auto strand: {Strand::FWD, Strand::REV}) {
      hashFilterCmpPair(mp2, p.complexity, chr.getStrandref(strand), p.complexity.getThreshold());
    }
    
    return;
  }
  
  /*int32_t check_sv(int32_t sv)
{
  // for paired-end
   LOG("   the read is paired in sequencing: %d\n",sv&1);
  LOG("   the read is mapped in a proper pair: %d\n",sv&2);
  LOG("   the query sequence itself is unmapped: %d\n",sv&4);
  LOG("   the mate is unmapped: %d\n",sv&8);
  LOG("   strand of the query (1 for reverse): %d\n",sv&16);
  LOG("   strand of the mate: %d\n",sv&32);
  LOG("   the read is the first read(F3) in a pair: %d\n",sv&64);
  LOG("   the read is the second read(F5) in a pair: %d\n",sv&128);
  LOG("   the alignment is not primary: %d\n",sv&256);
  LOG("   the read fails platform/vendor quality checks: %d\n",sv&512);
  LOG("   the read is either a PCR or an optical duplicate: %d\n",sv&1024);

  LOG("   template having multiple segments in sequencing: %d\n",sv&1);
  LOG("   each segment properly aligned according to the aligner: %d\n",sv&2);
  LOG("   segment unmapped: %d\n",sv&4);
  LOG("   next segment in the template unmapped: %d\n",sv&8);
  LOG("   SEQ being reverse complemented: %d\n",sv&16);
  LOG("   SEQ of the next segment in the template being reversed: %d\n",sv&32);
  LOG("   the first segment in the template: %d\n",sv&64);
  LOG("   the last segment in the template: %d\n",sv&128);
  LOG("   secondary alignment: %d\n",sv&256);
  LOG("   not passing quality controls: %d\n",sv&512);
  LOG("   PCR or optical duplicate: %d\n",sv&1024);
  LOG("   supplementary alignment: %d\n",sv&2048);


  // unmapped reads
  if(sv&4) goto err;
  // low quality reads
  if(sv&512 || sv&1024) goto err;
  //  if(p->rtype==READTYPE_PAIR){
    // unproper pair
    if(!(sv&2)) goto err;
    // unmatched pairs and interchromosomal pairs
    if(sv&8) goto err;
    // read pair mapped in same strand (for paired-end)
    if((sv&16 && sv&32) || (!(sv&16) && !(sv&32))) goto err;
    // }

 err:
  return 0;
}  */

}

void read_mapfile(const MyOpt::Variables &values, SeqStatsGenome &genome)
{
  std::vector<std::string> v;
  boost::split(v, values["input"].as<std::string>(), boost::algorithm::is_any_of(","));
  for(auto inputfile: v) {
    isFile(inputfile);
    std::cout << boost::format("Parsing %1%...\n") % inputfile;
    if(values.count("ftype")) {
      std::string ftype = values["ftype"].as<std::string>();
      if(ftype == "SAM" || ftype == "BAM") parseSam(values, inputfile, genome);
      else if(ftype == "BOWTIE") parseBowtie(values, inputfile, genome);
      else if(ftype == "TAGALIGN") parseTagAlign(values, inputfile, genome);
    } else {
      if(isStr(inputfile, ".sam") || isStr(inputfile, ".bam")) parseSam(values, inputfile, genome);
      else if(isStr(inputfile, ".bowtie"))   parseBowtie(values, inputfile, genome);
      else if(isStr(inputfile, ".tagalign")) parseTagAlign(values, inputfile, genome);
    }
  }

  if(!genome.getnread(Strand::BOTH)) PRINTERR("no read in input file.");

  return;
}

void checkRedundantReads(const MyOpt::Variables &values, Mapfile &p)
{
  p.complexity.setThreshold(values["thre_pb"].as<int32_t>(), p.genome);
  
  // Library complexity
  double r = getratio(values["ncmp"].as<int32_t>(), p.genome.getnread(Strand::BOTH));
  if(r>1){
    std::cerr << "Warning: number of reads is < "<< (int32_t)(values["ncmp"].as<int32_t>()/NUM_1M) <<" million.\n";
    p.complexity.lackOfRead_on();
  }
  p.complexity.setr4cmp(r*RAND_MAX);
  
  for(uint32_t i=0; i<p.genome.chr.size(); ++i) {
    if (values.count("pair")) filtering_eachchr_pair(p, p.genome.chr[i]);
    else                      filtering_eachchr_single(p, p.genome.chr[i]);
  }
  
  printf("done.\n");
  return;
}
