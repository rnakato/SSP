/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SHIFTPROFILE_H_
#define _SHIFTPROFILE_H_

#include <fstream>
#include <string>
#include <boost/bind.hpp>
#include "../common/inline.hpp"
#include "../common/BoostOptions.hpp"

class SeqStatsGenomeSSP;

class SSPstats {
  MyOpt::Opts opt;

  int32_t num4ssp;
  int32_t isExjac, isHd, isCc;
  int32_t eachchr;
  int32_t mp_from, mp_to;
  int32_t ng_from, ng_to, ng_step;
  int32_t numthreads;

  double nsc;
  double rlsc;
  double rsc;
  double backgroundUniformity;

 public:
  SSPstats() {}
  SSPstats(const int32_t ng_from_default,
	   const int32_t ng_to_default,
	   const int32_t ng_step_default,
	   const int32_t mp_from_default,
	   const int32_t mp_to_default):
  opt("Strand shift profile",100),
  num4ssp(0), isExjac(0), isHd(0), isCc(0), eachchr(0),
  mp_from(mp_from_default), mp_to(mp_to_default),
  ng_from(0), ng_to(0), ng_step(0), numthreads(0),
  nsc(0), rlsc(0), rsc(0), backgroundUniformity(0)
  {
    using namespace boost::program_options;
    opt.add_options()
      ("num4ssp",
       value<int32_t>()->default_value(NUM_10M)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, 1, "--num4ssp")),
       "Read number for calculating backgroud uniformity (per 100 Mbp)")
      ("ng_from",
       value<int32_t>()->default_value(ng_from_default)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, -1, "--ng_from")),
       "start shift of background")
      ("ng_to",
       value<int32_t>()->default_value(ng_to_default)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, -1, "--ng_to")),
       "end shift of background")
      ("ng_step",
       value<int32_t>()->default_value(ng_step_default)->notifier(boost::bind(&MyOpt::over<int32_t>, _1, -1, "--ng_step")),
       "step shift on of background")
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
    numthreads = MyOpt::getVal<int32_t>(values, "threads");

    DEBUGprint("SSPstats setValues done.");
  }
  void dump()
  {
    std::cout << boost::format("SSP background region: [%d,%d], step %d\n") % ng_from % ng_to % ng_step;
    std::cout << boost::format("Read number for SSP: %d\n") % num4ssp;
  }

  int32_t DoExjac()    const { return isExjac; }
  int32_t DoHd()       const { return isHd; }
  int32_t DoCc()       const { return isCc; }
  int32_t isEachchr()  const { return eachchr; }
  int32_t getnum4ssp() const { return num4ssp; }
  int32_t getMpFrom()  const { return mp_from; }
  int32_t getMpTo()    const { return mp_to; }
  int32_t getNgFrom()  const { return ng_from; }
  int32_t getNgTo()    const { return ng_to; }
  int32_t getNgStep()  const { return ng_step; }
  int32_t getnumthreads() const { return numthreads; }

  void setnsc(const double c) { nsc = c; }
  void setrlsc(const double c) { rlsc = c; }
  void setrsc(const double c) { rsc = c; }
  void setbu(const double c) { backgroundUniformity = c; }

  void printhead(std::ofstream &out) {
    out << "NSC\tRLSC\tRSC\tbackground uniformity";
  }
  void print(std::ofstream &out) {
    out << nsc << "\t" << rlsc << "\t" << rsc << "\t" << backgroundUniformity;
  }
};

template <class T, class S>
void makeProfile(SSPstats &sspst, S &genome, const std::string &head, const std::string &typestr);

void strShiftProfile(SSPstats &sspst, SeqStatsGenomeSSP &genome, const std::string &head, const std::string &typestr);

#endif /* _SHIFTPROFILE_H_ */
