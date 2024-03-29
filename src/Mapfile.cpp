/* Copyright(c) Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "Mapfile.hpp"

void FragmentLengthDist::outputDistFile(const std::string &prefix, const uint64_t nread) {
  std::string outputfile = prefix + ".ReadLengthDist.tsv";
  std::ofstream out(outputfile);
  printVector(out, vlenF3, "F3", nread);
  if(pairedend) {
    out << std::endl;
    printVector(out, vlenF5, "F5", nread);
  }
  out.close();

  if(pairedend) {
    outputfile = prefix + ".FragmentLengthDist.tsv";
    std::ofstream out(outputfile);
    printVector(out, vflen, "Fragmemt", nread);
  }
}

void SeqStatsGenomeSSP::setValues(const MyOpt::Variables &values) {
  DEBUGprint("SeqStatsGenomeSSP setValues...");

  inputfilename = MyOpt::getVal<std::string>(values, "input");
  pairedend     = values.count("pair");
  maxins        = MyOpt::getVal<int32_t>(values, "maxins");
  include_all_chromosome = values.count("include_allchr");
  specifyFtype  = values.count("ftype");
  if(onFtype()) {
    ftype = MyOpt::getVal<std::string>(values, "ftype");
    if(ftype != "SAM" && ftype != "BAM" && ftype != "CRAM" && ftype != "BOWTIE" && ftype != "TAGALIGN") PRINTERR_AND_EXIT("invalid --ftype.\n");
  }

  dflen.setValues(values);
  genometable = MyOpt::getVal<std::string>(values, "gt");
  readGenomeTable(genometable);

  if(values.count("mptable")) {
    for(auto &x: chr) x.getMptable(MyOpt::getVal<std::string>(values, "mptable"));
  }

  // Greekchr
  for(auto &x: chr) {
    if(x.getname() == "I") {
      for(auto &x:chr) x.Greekchron();
      break;
    }
  }

  // ignore isautosome();
  if(include_all_chromosome) {
    for(auto &x: chr) x.ConsiderAllchron();
  }
  
  // sepchr
  vsepchr = MyMthread::getVsepchr(getlen(), chr, MyOpt::getVal<int32_t>(values, "threads"));

  DEBUGprint("SeqStatsGenomeSSP setValues done.");
#ifdef DEBUG
  std::cout << "chr\tautosome" << std::endl;
  for(auto &x: chr) printList(x.getname(), x.isautosome());
  for(uint32_t i=0; i<vsepchr.size(); i++)
    std::cout << "thread " << (i+1) << ": " << vsepchr[i].s << "-" << vsepchr[i].e << std::endl;
  printReadstats();
#endif
}

void SeqStatsGenomeSSP::readGenomeTable(const std::string &gt)
{
  std::string lineStr;
  std::ifstream in(gt);
  if(!in) PRINTERR_AND_EXIT("Could not open " << gt << ".");

  printf("reading genome_table file..\n");

  try {
    while (!in.eof()) {
      std::vector<std::string> v;
      getline(in, lineStr);
      if(lineStr.empty() || lineStr[0] == '#') continue;
      if(ParseLine(v, lineStr, '\t')) PRINTERR_AND_EXIT("invalid format: " << gt );
      chr.emplace_back(v[0], stoi(v[1]));
//      std::cout << v[0] << ",,," << v[1] << std::endl;
    }
  } catch (const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    PRINTERR_AND_EXIT("invalid format of the genome_table file: " << gt );
  }

  return;
}
