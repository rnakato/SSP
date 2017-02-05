/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SSP_SHIFTPROFILE_H_
#define _SSP_SHIFTPROFILE_H_

#include <fstream>
#include <string>
#include <boost/bind.hpp>
#include "../common/inline.hpp"
#include "../common/BoostOptions.hpp"

class SeqStatsGenome;

class SSPstats {
  MyOpt::Opts opt;

  int32_t isExjac, isHd, isCc;
  int32_t eachchr;
  int32_t num4ssp;
  int32_t ng_from, ng_to, ng_step;
  int32_t ng_from_fcs, ng_to_fcs, ng_step_fcs;
  int32_t numthreads;
  
  double nsc;
  double rsc;
  double backgroundUniformity;
  double fcsread;
  double fcsflen;
  double fcs1k;
  double fcs10k;
  double fcs100k;

 public:
 SSPstats():
  opt("Strand shift profile",100),
    nsc(0), rsc(0), backgroundUniformity(0), fcsflen(0), fcs1k(0), fcs10k(0), fcs100k(0) {
    using namespace boost::program_options;
    opt.add_options()
      ("num4ssp",
       value<int32_t>()->default_value(NUM_10M)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--num4ssp")),
       "Read number for calculating backgroud uniformity (per 100 Mbp)")
      ("ng_from",
       value<int32_t>()->default_value(5*NUM_100K)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_from")),
       "start shift of background")
      ("ng_to",
       value<int32_t>()->default_value(NUM_1M)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_to")),
       "end shift of background")
      ("ng_step",
       value<int32_t>()->default_value(5000)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_step")),
       "step shift on of background")
      ("ng_from_fcs",
       value<int32_t>()->default_value(NUM_100K)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_from_fcs")),
       "fcs start of background")
      ("ng_to_fcs",
       value<int32_t>()->default_value(NUM_1M)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_to_fcs")),
       "fcs end of background")
      ("ng_step_fcs",
       value<int32_t>()->default_value(NUM_100K)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--ng_step_fcs")),
       "fcs step on of background")
      ("ssp_cc",    "make ssp based on cross correlation")
      ("ssp_hd",    "make ssp based on hamming distance")
      ("ssp_exjac", "make ssp based on extended Jaccard index")
      ("eachchr",   "make chromosome-sparated ssp files")
      ;
  }

  void setOpts(MyOpt::Opts &allopts) {
    allopts.add(opt);
  }
  
  void setValues(const MyOpt::Variables &values) {
    DEBUGprint("SSPstats setValues...");
    
    isExjac = values.count("ssp_exjac");
    isHd    = values.count("ssp_hd");
    isCc    = values.count("ssp_cc");
    eachchr = values.count("eachchr");
    num4ssp = MyOpt::getVal<int32_t>(values, "num4ssp");
    ng_from = MyOpt::getVal<int32_t>(values, "ng_from");
    ng_to   = MyOpt::getVal<int32_t>(values, "ng_to");
    ng_step = MyOpt::getVal<int32_t>(values, "ng_step");
    ng_from_fcs = MyOpt::getVal<int32_t>(values, "ng_from_fcs");
    ng_to_fcs   = MyOpt::getVal<int32_t>(values, "ng_to_fcs");
    ng_step_fcs = MyOpt::getVal<int32_t>(values, "ng_step_fcs");
    numthreads = MyOpt::getVal<int32_t>(values, "threads");
    
    DEBUGprint("SSPstats setValues done.");
  }
  void dump()
  {
    std::cout << boost::format("SSP background region: [%d,%d], step %d\n") % ng_from % ng_to % ng_step;
    std::cout << boost::format("FCS background region: [%d,%d], step %d\n") % ng_from_fcs % ng_to_fcs % ng_step_fcs;
    std::cout << boost::format("Read number for FCS: %d\n") % num4ssp;
  }

  int32_t DoExjac()    const { return isExjac; }
  int32_t DoHd()       const { return isHd; }
  int32_t DoCc()       const { return isCc; }
  int32_t isEachchr()  const { return eachchr; }
  int32_t getnum4ssp() const { return num4ssp; }
  int32_t getNgFrom()  const { return ng_from; }
  int32_t getNgTo()    const { return ng_to; }
  int32_t getNgStep()  const { return ng_step; }
  int32_t getNgFromFCS() const { return ng_from_fcs; }
  int32_t getNgToFCS()   const { return ng_to_fcs; }
  int32_t getNgStepFCS() const { return ng_step_fcs; }
  int32_t getnumthreads() const { return numthreads; }
  
  void setnsc(const double c) { nsc = c; }
  void setrsc(const double c) { rsc = c; }
  void setbu(const double c) { backgroundUniformity = c; }
  void setfcsread(const double c) { fcsread = c; }
  void setfcsflen(const double c) { fcsflen = c; }
  void setfcs1k(const double c) { fcs1k = c; }
  void setfcs10k(const double c) { fcs10k = c; }
  void setfcs100k(const double c) { fcs100k = c; }

  void printhead(std::ofstream &out) {
    out << "NSC\tRSC\tbackground uniformity\t"
	<< "FCS(read)\tFCS(flen)\tFCS(1k)\tFCS(10k)\tFCS(100k)"
	<< std::endl;
  }
  void print(std::ofstream &out) {
    out << nsc << "\t" << rsc << "\t" << backgroundUniformity << "\t" 
	<< fcsread << "\t" << fcsflen << "\t" << fcs1k << "\t" << fcs10k
	<< "\t" << fcs100k << std::endl;
  }
};


void strShiftProfile(SSPstats &sspst, SeqStatsGenome &genome, const std::string &head, const std::string &typestr);
void makeFCSProfile(SSPstats &sspst, const SeqStatsGenome &genome, const std::string &head, const std::string &typestr);

#endif /* _SSP_SHIFTPROFILE_H_ */
