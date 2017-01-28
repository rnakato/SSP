/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
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
#include "ssp_gv.h"
#include "pw_readmapfile.h"
#include "ssp_shiftprofile.h"

#define VERSION "1.0.0"

namespace {
  const int32_t numGcov(5000000);
}

MyOpt::Variables getOpts(SSP::Global &ssp, int argc, char* argv[]);
void setOpts(MyOpt::Opts &);
void init_dump(const MyOpt::Variables &);

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

void estimateFragLength(const MyOpt::Variables &values, SSP::Global &p)
{
  if(values.count("pair") || values.count("nomodel")) return;

  std::string head(p.getprefix());
  
  clock_t t1,t2;
  t1 = clock();
  strShiftProfile(p.sspst, values, p.genome, head, "jaccard"); 
  t2 = clock();
  std::cout << "Jaccard Bit: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  makeFCSProfile(p.sspst, values, p.genome, head, "fcs");
  p.outputSSPstats();

  clock_t t3 = clock();
  std::cout << "Fragment variability: " << static_cast<double>(t3 - t2) / CLOCKS_PER_SEC << "sec.\n";

  if(values.count("ssp_exjac")) {
    t1 = clock();
    strShiftProfile(p.sspst, values, p.genome, head, "exjaccard");
    t2 = clock();
    std::cout << "Jaccard Vec: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }

  if(values.count("ssp_hd")) {
    t1 = clock();
    strShiftProfile(p.sspst, values, p.genome, head, "hdp");
    t2 = clock();
    std::cout << "Hamming: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }
    
  if(values.count("ssp_cc")) {
    t1 = clock();
    strShiftProfile(p.sspst, values, p.genome, head, "ccp");
    t2 = clock();    
    std::cout << "ccp: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }
  
  return;
}

int main(int argc, char* argv[])
{
  SSP::Global p;
  MyOpt::Variables values = getOpts(p, argc, argv);
  p.setValues(values);

  read_mapfile(values, p.genome);
  p.complexity.checkRedundantReads(p.genome);
  
  estimateFragLength(values, p);

#ifdef DEBUG
  p.genome.printReadstats();
#endif

  return 0;
}

void checkParam(const MyOpt::Variables &values)
{
  DEBUGprint("checkParam...");

  if(values.count("ftype")) {
    std::string ftype = values["ftype"].as<std::string>();
    if(ftype != "SAM" && ftype != "BAM" && ftype != "BOWTIE" && ftype != "TAGALIGN") PRINTERR("invalid --ftype.\n");
  }
  
  DEBUGprint("checkParam done.");
  return;
}

MyOpt::Variables getOpts(SSP::Global &ssp, int argc, char* argv[])
{
  DEBUGprint("setOpts...");

  MyOpt::Opts allopts("Options");

  ssp.setOpts(allopts);
  
  setOpts(allopts);
  
  MyOpt::Variables values;
  
  DEBUGprint("getOpts...");

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
      std::cout << "\n" << allopts << std::endl;
      exit(0);
    }
    std::vector<std::string> opts = {"input", "output", "gt"};
    for (auto x: opts) {
      if (!values.count(x)) PRINTERR("specify --" << x << " option.");
    }

    notify(values);
    checkParam(values);
    
    boost::filesystem::path dir(values["odir"].as<std::string>());
    boost::filesystem::create_directory(dir);
  
    init_dump(values);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  
  DEBUGprint("getOpts done.");
  return values;
}

void setOpts(MyOpt::Opts &allopts)
{
  using namespace boost::program_options;

  MyOpt::Opts optIO("Input/Output",100);
  optIO.add_options()
    ("input,i",   value<std::string>(), "Mapping file. Multiple files are allowed (separated by ',')")
    ("output,o",  value<std::string>(), "Prefix of output files")
    ("odir",    value<std::string>()->default_value("sspout"), "output directory name")
    ("pair", "add when the input file is paired-end")
    ;
  MyOpt::Opts optOptional("Optional",100);
  optOptional.add_options()
    ("maxins",     value<int32_t>()->default_value(500), "maximum fragment length")
    ("ftype,f", value<std::string>(), "{SAM|BAM|BOWTIE|TAGALIGN}: format of input file\nTAGALIGN could be gzip'ed (extension: tagAlign.gz)")
    ;
  MyOpt::Opts optOther("Others",100);
  optOther.add_options()
    ("threads,p", value<int32_t>()->default_value(1)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--thread")),
     "number of threads to launch")
    ("version,v", "print version")
    ("help,h", "show help message")
    ;
  allopts.add(optIO).add(optOptional).add(optOther);
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
