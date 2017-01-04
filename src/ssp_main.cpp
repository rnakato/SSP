/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "pw_readmapfile.h"
#include "ssp_estFlen.h"
#include "pw_gv.h"

#define VERSION "1.0.0"

namespace {
  const int numGcov(5000000);
}

MyOpt::Variables getOpts(int argc, char* argv[]);
void setOpts(MyOpt::Opts &, MyOpt::Opts &);
void init_dump(const MyOpt::Variables &);
void output_stats(const MyOpt::Variables &values, const Mapfile &p);
void output_wigstats(Mapfile &p);

void SeqStatsGenome::readGenomeTable(const std::string &gt, const int binsize)
{
  std::vector<std::string> v;
  std::string lineStr;
  std::ifstream in(gt);
  if(!in) PRINTERR("Could nome open " << gt << ".");

  while (!in.eof()) {
    getline(in, lineStr);
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    SeqStats s(v[0], stoi(v[1]));
    s.nbin = s.getlen()/binsize +1;
    chr.push_back(s);
  }

  return;
}

std::vector<sepchr> SeqStatsGenome::getVsepchr(const int numthreads)
{
  std::vector<sepchr> vsepchr;

  uint sepsize = len/numthreads;
  for(uint i=0; i<chr.size(); ++i) {
    uint s = i;
    long len(0);
    while(len < sepsize && i<chr.size()) {
      len += chr[i].getlen();
      i++;
    }
    i--;
    uint e = i;
    sepchr sep(s,e);
    vsepchr.push_back(sep);
  }
  return vsepchr;
}

void printVersion()
{
  std::cerr << "SSP version " << VERSION << std::endl;
  exit(0);
}

void help_global()
{
  auto helpmsg = R"(
===============

Usage: ssp [option] -i <inputfile> -o <output> --gt <genome_table>)";
  
  std::cerr << "\nSSP v" << VERSION << helpmsg << std::endl;
  return;
}

int main(int argc, char* argv[])
{
  MyOpt::Variables values = getOpts(argc, argv);

  boost::filesystem::path dir(values["odir"].as<std::string>());
  boost::filesystem::create_directory(dir);

  Mapfile p(values);
  read_mapfile(values, p);

  if(!values.count("nofilter")) checkRedundantReads(values, p);
  else                          p.genome.setnread2nread_red();
  p.genome.setnread_red();
  
  estimateFragLength(values, p);

#ifdef DEBUG
  p.printstats();
#endif

  //  output_stats(values, p);

  return 0;
}

void checkParam(const MyOpt::Variables &values)
{
#ifdef DEBUG
  std::cout << "checkParam..." << std::endl;
#endif

  std::vector<std::string> intopts = {"threads"};
  for (auto x: intopts) chkminus<int>(values, x, 0);
  std::vector<std::string> intopts2 = {"thre_pb"};
  for (auto x: intopts2) chkminus<int>(values, x, -1);

  if(values.count("ftype")) {
    std::string ftype = values["ftype"].as<std::string>();
    if(ftype != "SAM" && ftype != "BAM" && ftype != "BOWTIE" && ftype != "TAGALIGN") PRINTERR("invalid --ftype.\n");
  }
  
#ifdef DEBUG
  std::cout << "checkParam done." << std::endl;
#endif
  return;
}

MyOpt::Variables getOpts(int argc, char* argv[])
{
#ifdef DEBUG
  std::cout << "setOpts..." << std::endl;
#endif

  MyOpt::Opts allopts("Options");
  MyOpt::Opts opts4help("Options");
  setOpts(allopts, opts4help);
  
  MyOpt::Variables values;
  
#ifdef DEBUG
  std::cout << "getOpts..." << std::endl;
#endif

  try {
    boost::program_options::parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);
    if (values.count("version")) printVersion();
    if (argc ==1) {
      help_global();
      std::cerr << "Use --help option for more information on the other options\n\n";
      exit(0);
    }
    if (values.count("help")) {
      help_global();
      std::cout << "\n" << opts4help << std::endl;
      exit(0);
    }
    std::vector<std::string> opts = {"input", "output", "gt"};
    for (auto x: opts) {
      if (!values.count(x)) PRINTERR("specify --" << x << " option.");
    }

    notify(values);
    checkParam(values);
  
    init_dump(values);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  
#ifdef DEBUG
  std::cout << "getOpts done." << std::endl;
#endif
  return values;
}

void setOpts(MyOpt::Opts &allopts,MyOpt::Opts &opts4help)
{
  using namespace boost::program_options;

  MyOpt::Opts optreq("Required",100);
  optreq.add_options()
    ("input,i",   value<std::string>(), "Mapping file. Multiple files are allowed (separated by ',')")
    ("output,o",  value<std::string>(), "Prefix of output files")
    ("gt",        value<std::string>(), "Genome table (tab-delimited file describing the name and length of each chromosome)")
    ;
  
  MyOpt::Opts optssp("Strand shift profile",100);
  optssp.add_options()
    ("ng_from", value<int>()->default_value(5*NUM_100K), "start shift of background")
    ("ng_to",   value<int>()->default_value(NUM_1M),     "end shift of background")
    ("ng_step", value<int>()->default_value(5000),       "step shift on of background")
    ("num4ssp", value<int>()->default_value(NUM_10M),    "Read number for calculating backgroud uniformity (per 100 Mbp)")
    ("ssp_cc",    "make ssp based on cross correlation")
    ("ssp_hd",    "make ssp based on hamming distance")
    ("ssp_exjac", "make ssp based on extended Jaccard index")
    ("eachchr", "make chromosome-sparated ssp files")
    ("mptable", value<std::string>(), "Genome table for mappable regions")
    ;
  MyOpt::Opts optfcs("Fragment cluster score",100);
  optfcs.add_options()
    ("ng_from_fcs", value<int>()->default_value(NUM_100K), "fcs start of background")
    ("ng_to_fcs",   value<int>()->default_value(NUM_1M),   "fcs end of background")
    ("ng_step_fcs", value<int>()->default_value(NUM_100K), "fcs step on of background")
    ;
  MyOpt::Opts optIO("Optional",100);
  optIO.add_options()
    ("ftype,f", value<std::string>(), "{SAM|BAM|BOWTIE|TAGALIGN}: format of input file\nTAGALIGN could be gzip'ed (extension: tagAlign.gz)")
    ("odir",    value<std::string>()->default_value("sspout"),	  "output directory name")
    ("fcsfull",  "outout full fragment variability profile")
    ("nofilter", "do not filter PCR bias")
    ("thre_pb", value<int>()->default_value(0),	       "PCRbias threshold (default: more than max(1 read, 10 times greater than genome average)) ")
    ("ncmp",    value<int>()->default_value(10000000), "read number for calculating library complexity")
    ;
  MyOpt::Opts optother("Others",100);
  optother.add_options()
    ("threads,p",    value<int>()->default_value(1),  "number of threads to launch")
    ("version,v", "print version")
    ("help,h", "show help message")
    ;

  // Ignored options
  MyOpt::Opts optignore("for parse2wig",100);
  optignore.add_options()
    ("mp",        value<std::string>(),	  "Mappability file")
    ("mpthre",    value<double>()->default_value(0.3),	  "Threshold of low mappability regions")
    ("flen",        value<int>()->default_value(150), "predefined fragment length\n(Automatically calculated in paired-end mode)")
    ("nomodel",   "predefine the fragment length (default: estimated by hamming distance plot)")
    ("binsize,b",   value<int>()->default_value(50),	  "bin size")
    ("of",        value<int>()->default_value(0),	  "output format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
    ("rcenter", value<int>()->default_value(0), "consider length around the center of fragment ")
    ("pair", 	  "add when the input file is paired-end")
    ("maxins",     value<int>()->default_value(500), "maximum fragment length")
    ("genome",     value<std::string>(),	  "reference genome sequence for GC content estimation")
    ("flen4gc",    value<int>()->default_value(120),  "fragment length for calculation of GC distribution")
    ("gcdepthoff", "do not consider depth of GC contents")
    ("ntype,n",        value<std::string>()->default_value("NONE"),  "Total read normalization\n{NONE|GR|GD|CR|CD}\n   NONE: not normalize\n   GR: for whole genome, read number\n   GD: for whole genome, read depth\n   CR: for each chromosome, read number\n   CD: for each chromosome, read depth")
    ("nrpm",        value<int>()->default_value(20000000),	  "Total read number after normalization")
    ("ndepth",      value<double>()->default_value(1.0),	  "Averaged read depth after normalization")
    ("bed",        value<std::string>(),	  "specify the BED file of enriched regions (e.g., peak regions)")
    ;  
    ;
    
  allopts.add(optreq).add(optssp).add(optfcs).add(optIO).add(optother).add(optignore);
  opts4help.add(optreq).add(optssp).add(optfcs).add(optIO).add(optother);
  return;
}

void init_dump(const MyOpt::Variables &values){
  std::vector<std::string> str_wigfiletype = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
 
  BPRINT("\n======================================\n");
  BPRINT("SSP version %1%\n\n") % VERSION;
  BPRINT("Input file %1%\n")         % values["input"].as<std::string>();
  if(values.count("ftype")) BPRINT("\tFormat: %1%\n") % values["ftype"].as<std::string>();
  BPRINT("Output file: %1%/%2%\n")   % values["odir"].as<std::string>() % values["output"].as<std::string>();
  BPRINT("Genome-table file: %1%\n") % values["gt"].as<std::string>();
  if(values.count("mptable")) BPRINT("Mappable genome-table file: %1%\n") % values["mptable"].as<std::string>();
  
  if (!values.count("nofilter")) {
    BPRINT("PCR bias filtering: ON\n");
    if (values["thre_pb"].as<int>()) BPRINT("PCR bias threshold: > %1%\n") % values["thre_pb"].as<int>();
  } else {
    BPRINT("PCR bias filtering: OFF\n");
  }
  if(values.count("fcsfull")) BPRINT("\tplot full fcs profile.\n");
  BPRINT("num4ssp %d\n") % values["num4ssp"].as<int>();
  BPRINT("background region: [%d,%d], step %d\n") % values["ng_from"].as<int>() % values["ng_to"].as<int>() % values["ng_step"].as<int>();
  
  BPRINT("\nNumber of threads: %1%\n") % values["threads"].as<int>();
  printf("======================================\n");
  return;
}

void print_SeqStats(const MyOpt::Variables &values, std::ofstream &out, const SeqStats &p, const Mapfile &mapfile)
{
  /* genome data */
  out << p.name << "\t" << p.getlen()  << "\t" << p.getlenmpbl() << "\t" << p.getpmpbl() << "\t";
  /* total reads*/
  out << boost::format("%1%\t%2%\t%3%\t%4$.1f%%\t")
    % p.bothnread() % p.seq[STRAND_PLUS].nread % p.seq[STRAND_MINUS].nread
    % (p.bothnread()*100/static_cast<double>(mapfile.genome.bothnread()));

  /* nonredundant reads */
  printr(out, p.bothnread_nonred(), p.bothnread());
  p.seq[STRAND_PLUS].printnonred(out);
  p.seq[STRAND_MINUS].printnonred(out);
  printr(out, p.bothnread_red(), p.bothnread());
  p.seq[STRAND_PLUS].printred(out);
  p.seq[STRAND_MINUS].printred(out);

  /* reads after GCnorm */
  if(values.count("genome")) {
    printr(out, p.bothnread_afterGC(), p.bothnread());
    p.seq[STRAND_PLUS].printafterGC(out);
    p.seq[STRAND_MINUS].printafterGC(out);
  }
  out << boost::format("%1$.3f\t") % p.depth;
  if(p.getweight4rpm()) out << boost::format("%1$.3f\t") % p.getweight4rpm(); else out << " - \t";
  if(values["ntype"].as<std::string>() == "NONE") out << p.bothnread_nonred() << "\t"; else out << p.bothnread_rpm() << "\t";

  p.printGcov(out, mapfile.islackOfRead4GenomeCov());
  
  p.ws.printPoispar(out);
  if(values.count("bed")) out << boost::format("%1$.3f\t") % p.getFRiP();

  p.ws.printZINBpar(out);
  
  out << std::endl;
  return;
}

void output_stats(const MyOpt::Variables &values, const Mapfile &p)
{
  std::string filename = p.getbinprefix() + ".csv";
  std::ofstream out(filename);

  out << "SSP version " << VERSION << std::endl;
  out << "Input file: \"" << values["input"].as<std::string>() << "\"" << std::endl;
  out << "Redundancy threshold: >" << p.getthre4filtering() << std::endl;

  p.printComplexity(out);
  p.printFlen(values, out);
  if(values.count("genome")) out << "GC summit: " << p.getmaxGC() << std::endl;

  // Global stats
  out << "\n\tlength\tmappable base\tmappability\t";
  out << "total reads\t\t\t\t";
  out << "nonredundant reads\t\t\t";
  out << "redundant reads\t\t\t";
  if(values.count("genome")) out << "reads (GCnormed)\t\t\t";
  out << "read depth\t";
  out << "scaling weight\t";
  out << "normalized read number\t";
  out << "gcov (Raw)\tgcov (Normed)\t";
  out << "bin mean\tbin variance\t";
  if(values.count("bed")) out << "FRiP\t";
  out << "nb_p\tnb_n\tnb_p0\t";
  out << std::endl;
  out << "\t\t\t\t";
  out << "both\tforward\treverse\t% genome\t";
  out << "both\tforward\treverse\t";
  out << "both\tforward\treverse\t";
  if(values.count("genome")) out << "both\tforward\treverse\t";
  out << std::endl;

  // SeqStats
  print_SeqStats(values, out, p.genome, p);
  for(auto x:p.genome.chr) print_SeqStats(values, out, x, p);
  
  std::cout << "stats is output in " << filename << "." << std::endl;

  return;
}

std::vector<char> makeGcovArray(const MyOpt::Variables &values, SeqStats &chr, Mapfile &p, double r4cmp)
{
  std::vector<char> array;
  if(values.count("mp")) array = readMpbl_binary(values["mp"].as<std::string>(), ("chr" + p.lchr->name), chr.getlen());
  else array = readMpbl_binary(chr.getlen());
  if(values.count("bed")) arraySetBed(array, chr.name, p.genome.getvbed());

  int val(0);
  int size = array.size();
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto &x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      
      if(rand() >= r4cmp) val=COVREAD_ALL; else val=COVREAD_NORM;
      
      int s(std::max(0, std::min(x.F3, x.F5)));
      int e(std::min(std::max(x.F3, x.F5), size-1));
      if(s >= size || e < 0) {
	std::cerr << "Warning: " << chr.name << " read " << s <<"-"<< e << " > array size " << array.size() << std::endl;
      }
      for(int i=s; i<=e; ++i) if(array[i]==MAPPABLE) array[i]=val;
    }
  }
  return array;
}

void calcGcovchr(const MyOpt::Variables &values, Mapfile &p, int s, int e, double r4cmp, boost::mutex &mtx)
{
  for(int i=s; i<=e; ++i) {
    std::cout << p.genome.chr[i].name << ".." << std::flush;
    auto array = makeGcovArray(values, p.genome.chr[i], p, r4cmp);
    p.genome.chr[i].calcGcov(array);
    p.genome.addGcov(i, mtx);
  }
}

void calcFRiP(SeqStats &chr, const std::vector<bed> vbed)
{
  std::vector<char> array(chr.getlen(), MAPPABLE);
  arraySetBed(array, chr.name, vbed);
  for(int strand=0; strand<STRANDNUM; ++strand) {
    for (auto &x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      int s(std::min(x.F3, x.F5));
      int e(std::max(x.F3, x.F5));
      for(int i=s; i<=e; ++i) {
	if(array[i]==INBED) {
	  x.inpeak = 1;
	  ++chr.nread_inbed;
	  break;
	}
      }
    }
  }
  return;
}

void output_wigstats(Mapfile &p)
{
  std::string filename = p.getbinprefix() + ".binarray_dist.csv";
  std::ofstream out(filename);

  std::cout << "generate " << filename << ".." << std::flush;

  out << "\tGenome\t\t\t";
  for (auto x:p.genome.chr) out << x.name << "\t\t\t\t";
  out << std::endl;
  out << "read number\tnum of bins genome\tprop\tZINB estimated\t";
  for (auto x:p.genome.chr) out << "num of bins\tprop\tPoisson estimated\tZINB estimated\t";
  out << std::endl;

  for(size_t i=0; i<p.genome.ws.wigDist.size(); ++i) {
    out << i << "\t";
    p.genome.ws.printwigDist(out, i);
    out << p.genome.ws.getZINB(i) << "\t";
    //    out << p.genome.getZIP(i) << "\t";
    for (auto x:p.genome.chr) {
      x.ws.printwigDist(out, i);
      out << x.ws.getPoisson(i) << "\t";
      out << x.ws.getZINB(i) << "\t";
    }
    out << std::endl;
  }

  std::cout << "done." << std::endl;
  return;
}
