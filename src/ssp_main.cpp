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
  read_mapfile(values, p.genome);

  if(!values.count("nofilter")) checkRedundantReads(values, p);
  
  estimateFragLength(values, p);

#ifdef DEBUG
  p.genome.printReadstats();
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
 
  std::cout << boost::format("\n======================================\n");
  std::cout << boost::format("SSP version %1%\n\n") % VERSION;
  std::cout << boost::format("Input file %1%\n")         % values["input"].as<std::string>();
  if(values.count("ftype")) std::cout << boost::format("\tFormat: %1%\n") % values["ftype"].as<std::string>();
  std::cout << boost::format("Output file: %1%/%2%\n")   % values["odir"].as<std::string>() % values["output"].as<std::string>();
  std::cout << boost::format("Genome-table file: %1%\n") % values["gt"].as<std::string>();
  if(values.count("mptable")) std::cout << boost::format("Mappable genome-table file: %1%\n") % values["mptable"].as<std::string>();
  
  if (!values.count("nofilter")) {
    std::cout << boost::format("PCR bias filtering: ON\n");
    if (values["thre_pb"].as<int32_t>()) std::cout << boost::format("PCR bias threshold: > %1%\n") % values["thre_pb"].as<int32_t>();
  } else {
    std::cout << boost::format("PCR bias filtering: OFF\n");
  }
  std::cout << boost::format("num4ssp %d\n") % values["num4ssp"].as<int32_t>();
  std::cout << boost::format("background region: [%d,%d], step %d\n")
    % values["ng_from"].as<int32_t>()
    % values["ng_to"].as<int32_t>()
    % values["ng_step"].as<int32_t>();
  
  std::cout << boost::format("\nNumber of threads: %1%\n") % values["threads"].as<int32_t>();
  printf("======================================\n");
  return;
}

void estimateFragLength(const MyOpt::Variables &values, Mapfile &p)
{
  if(values.count("pair") || values.count("nomodel")) return;

  clock_t t1,t2;
  t1 = clock();
  strShiftProfile(values, p, "jaccard"); 
  t2 = clock();
  std::cout << "Jaccard Bit: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  makeFCSProfile(values, p, "fcs");
  p.outputSSPstats();

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
