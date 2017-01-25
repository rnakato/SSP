/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "pw_gv.h"
#include "pw_readmapfile.h"
#include "ssp_shiftprofile.h"

#define VERSION "1.0.0"

namespace {
  const int32_t numGcov(5000000);
}

MyOpt::Variables getOpts(int argc, char* argv[]);
void setOpts(MyOpt::Opts &, MyOpt::Opts &);
void init_dump(const MyOpt::Variables &);
void output_stats(const MyOpt::Variables &values, const Mapfile &p);
//void output_wigstats(Mapfile &p);
void estimateFragLength(const MyOpt::Variables &values, Mapfile &p);

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
  //  else                          p.genome.setnread2nread_red();
  //  p.genome.setnread_red();
  
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
  for (auto x: intopts) chkminus<int32_t>(values, x, 0);
  std::vector<std::string> intopts2 = {"thre_pb"};
  for (auto x: intopts2) chkminus<int32_t>(values, x, -1);

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
    ("ng_from", value<int32_t>()->default_value(5*NUM_100K), "start shift of background")
    ("ng_to",   value<int32_t>()->default_value(NUM_1M),     "end shift of background")
    ("ng_step", value<int32_t>()->default_value(5000),       "step shift on of background")
    ("num4ssp", value<int32_t>()->default_value(NUM_10M),    "Read number for calculating backgroud uniformity (per 100 Mbp)")
    ("ssp_cc",    "make ssp based on cross correlation")
    ("ssp_hd",    "make ssp based on hamming distance")
    ("ssp_exjac", "make ssp based on extended Jaccard index")
    ("eachchr", "make chromosome-sparated ssp files")
    ("mptable", value<std::string>(), "Genome table for mappable regions")
    ;
  MyOpt::Opts optfcs("Fragment cluster score",100);
  optfcs.add_options()
    ("ng_from_fcs", value<int32_t>()->default_value(NUM_100K), "fcs start of background")
    ("ng_to_fcs",   value<int32_t>()->default_value(NUM_1M),   "fcs end of background")
    ("ng_step_fcs", value<int32_t>()->default_value(NUM_100K), "fcs step on of background")
    ;
  MyOpt::Opts optIO("Optional",100);
  optIO.add_options()
    ("ftype,f", value<std::string>(), "{SAM|BAM|BOWTIE|TAGALIGN}: format of input file\nTAGALIGN could be gzip'ed (extension: tagAlign.gz)")
    ("odir",    value<std::string>()->default_value("sspout"),	  "output directory name")
    ("nofilter", "do not filter PCR bias")
    ("thre_pb", value<int32_t>()->default_value(0),	       "PCRbias threshold (default: more than max(1 read, 10 times greater than genome average)) ")
    ("ncmp",    value<int32_t>()->default_value(10000000), "read number for calculating library complexity")
    ;
  MyOpt::Opts optother("Others",100);
  optother.add_options()
    ("threads,p",    value<int32_t>()->default_value(1),  "number of threads to launch")
    ("version,v", "print version")
    ("help,h", "show help message")
    ;

  // Ignored options
  MyOpt::Opts optignore("for parse2wig",100);
  optignore.add_options()
    ("mp",        value<std::string>(),	  "Mappability file")
    ("mpthre",    value<double>()->default_value(0.3),	  "Threshold of low mappability regions")
    ("flen",        value<int32_t>()->default_value(150), "predefined fragment length\n(Automatically calculated in paired-end mode)")
    ("nomodel",   "predefine the fragment length (default: estimated by hamming distance plot)")
    ("binsize,b",   value<int32_t>()->default_value(50),	  "bin size")
    ("of",        value<int32_t>()->default_value(0),	  "output format\n   0: binary (.bin)\n   1: compressed wig (.wig.gz)\n   2: uncompressed wig (.wig)\n   3: bedGraph (.bedGraph)\n   4: bigWig (.bw)")
    ("rcenter", value<int32_t>()->default_value(0), "consider length around the center of fragment ")
    ("pair", 	  "add when the input file is paired-end")
    ("maxins",     value<int32_t>()->default_value(500), "maximum fragment length")
    ("genome",     value<std::string>(),	  "reference genome sequence for GC content estimation")
    ("flen4gc",    value<int32_t>()->default_value(120),  "fragment length for calculation of GC distribution")
    ("gcdepthoff", "do not consider depth of GC contents")
    ("ntype,n",        value<std::string>()->default_value("NONE"),  "Total read normalization\n{NONE|GR|GD|CR|CD}\n   NONE: not normalize\n   GR: for whole genome, read number\n   GD: for whole genome, read depth\n   CR: for each chromosome, read number\n   CD: for each chromosome, read depth")
    ("nrpm",        value<int32_t>()->default_value(20000000),	  "Total read number after normalization")
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
    if (values["thre_pb"].as<int32_t>()) BPRINT("PCR bias threshold: > %1%\n") % values["thre_pb"].as<int32_t>();
  } else {
    BPRINT("PCR bias filtering: OFF\n");
  }
  BPRINT("num4ssp %d\n") % values["num4ssp"].as<int32_t>();
  BPRINT("background region: [%d,%d], step %d\n")
    % values["ng_from"].as<int32_t>()
    % values["ng_to"].as<int32_t>()
    % values["ng_step"].as<int32_t>();
  
  BPRINT("\nNumber of threads: %1%\n") % values["threads"].as<int32_t>();
  printf("======================================\n");
  return;
}

void print_SeqStats(const MyOpt::Variables &values, std::ofstream &out, const SeqStats &p, const Mapfile &mapfile)
{
  /* genome data */
  out << p.name << "\t" << p.getlen()  << "\t" << p.getlenmpbl() << "\t" << p.getpmpbl() << "\t";
  /* total reads*/
  out << boost::format("%1%\t%2%\t%3%\t%4$.1f%%\t")
    % p.getnread(STRAND_BOTH) % p.getnread(STRAND_PLUS) % p.getnread(STRAND_MINUS)
    % (p.getnread(STRAND_BOTH)*100/static_cast<double>(mapfile.genome.getnread(STRAND_BOTH)));

  /* nonredundant reads */
  printr(out, p.getnread_nonred(STRAND_BOTH), p.getnread(STRAND_BOTH));
  //  p.seq[STRAND_PLUS].printnonred(out);
  // p.seq[STRAND_MINUS].printnonred(out);
  printr(out, p.getnread_red(STRAND_BOTH), p.getnread(STRAND_BOTH));
  //  p.seq[STRAND_PLUS].printred(out);
  // p.seq[STRAND_MINUS].printred(out);

  /* reads after GCnorm */
  if(values.count("genome")) {
    printr(out, p.getnread_afterGC(STRAND_BOTH), p.getnread(STRAND_BOTH));
    //  p.seq[STRAND_PLUS].printafterGC(out);
    // p.seq[STRAND_MINUS].printafterGC(out);
  }
  out << boost::format("%1$.3f\t") % p.getdepth();
  if(p.getweight4rpm()) out << boost::format("%1$.3f\t") % p.getweight4rpm(); else out << " - \t";
  if(values["ntype"].as<std::string>() == "NONE") out << p.getnread_nonred(STRAND_BOTH) << "\t";
  else out << p.getnread_rpm(STRAND_BOTH) << "\t";

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
  for(auto &x:p.genome.chr) print_SeqStats(values, out, x, p);
  
  std::cout << "stats is output in " << filename << "." << std::endl;

  return;
}

void calcFRiP(SeqStats &chr, const std::vector<bed> vbed)
{
  std::vector<int8_t> array(chr.getlen(), MAPPABLE);
  arraySetBed(array, chr.name, vbed);
  for(int32_t strand=0; strand<STRANDNUM; ++strand) {
    for (auto &x: chr.seq[strand].vRead) {
      if(x.duplicate) continue;
      int32_t s(std::min(x.F3, x.F5));
      int32_t e(std::max(x.F3, x.F5));
      for(int32_t i=s; i<=e; ++i) {
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

/*void output_wigstats(Mapfile &p)
{
  std::string filename = p.getbinprefix() + ".binarray_dist.csv";
  std::ofstream out(filename);

  std::cout << "generate " << filename << ".." << std::flush;

  out << "\tGenome\t\t\t";
  for (auto &x: p.genome.chr) out << x.name << "\t\t\t\t";
  out << std::endl;
  out << "read number\tnum of bins genome\tprop\tZINB estimated\t";
  for (auto &x: p.genome.chr) out << "num of bins\tprop\tPoisson estimated\tZINB estimated\t";
  out << std::endl;

  for(size_t i=0; i<p.genome.ws.wigDist.size(); ++i) {
    out << i << "\t";
    p.genome.ws.printwigDist(out, i);
    out << p.genome.ws.getZINB(i) << "\t";
    //    out << p.genome.getZIP(i) << "\t";
    for (auto &x:p.genome.chr) {
      x.ws.printwigDist(out, i);
      out << x.ws.getPoisson(i) << "\t";
      out << x.ws.getZINB(i) << "\t";
    }
    out << std::endl;
  }

  std::cout << "done." << std::endl;
  return;
  }*/

void estimateFragLength(const MyOpt::Variables &values, Mapfile &p)
{
  if(values.count("pair") || values.count("nomodel")) return;

  clock_t t1,t2;
  t1 = clock();
  strShiftProfile(values, p, "jaccard"); 
  t2 = clock();
  std::cout << "Jaccard Bit: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  makeFCSProfile(values, p, "fcs");
  p.printSSPstats();

  clock_t t3 = clock();
  std::cout << "Fragment variability: " << static_cast<double>(t3 - t2) / CLOCKS_PER_SEC << "sec.\n";

  if(values.count("ssp_exjac")) {
    t1 = clock();
    strShiftProfile(values, p, "exjaccard");
    t2 = clock();
    std::cout << "Jaccard Vec: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }

  if(values.count("ssp_hd")) {
    t1 = clock();
    strShiftProfile(values, p, "hdp");
    t2 = clock();
    std::cout << "Hamming: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }
    
  if(values.count("ssp_cc")) {
    t1 = clock();
    strShiftProfile(values, p, "ccp");
    t2 = clock();    
    std::cout << "ccp: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }
  
  return;
}
