/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <boost/bind.hpp>
#include "BoostOptions.hpp"

namespace MyOpt {
  void setOptIO(Opts &allopts, const std::string &odir_default)
  {
    Opts opt("Input/Output",100);
    opt.add_options()
      ("input,i", boost::program_options::value<std::string>(),
       "Mapping file. Multiple files are allowed (separated by ',')")
      ("output,o",boost::program_options::value<std::string>(),
       "Prefix of output files")
      ("odir",    boost::program_options::value<std::string>()->default_value(odir_default),
       "output directory name")
      ("ftype,f", boost::program_options::value<std::string>(),
       "{SAM|BAM|BOWTIE|TAGALIGN}: format of input file\nTAGALIGN could be gzip'ed (extension: tagAlign.gz)")
      ;
    allopts.add(opt);
  }
  void setOptPair(Opts &allopts)
  {
    MyOpt::Opts opt("For paired-end",100);
    opt.add_options()
      ("pair",   "add when the input file is paired-end")
      ("maxins",
       boost::program_options::value<int32_t>()->default_value(500)->notifier(boost::bind(&over<int32_t>, _1, 1, "--maxins")),
       "maximum fragment length")
      ;
    allopts.add(opt);
  }
  void setOptOther(Opts &allopts)
  {
    Opts opt("Others",100);
    opt.add_options()
      ("threads,p",
       boost::program_options::value<int32_t>()->default_value(1)->notifier(boost::bind(&over<int32_t>, _1, 1, "--thread")),
       "number of threads to launch")
      ("version,v", "print version")
      ("help,h", "show help message")
      ;
    allopts.add(opt);
  }
  
}
