/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SSP_GV_H_
#define _SSP_GV_H_

#include <fstream>
#include <numeric>
#include <boost/thread.hpp>
#include "readdata.h"
#include "mthread.h"
#include "LibraryComplexity.hpp"
#include "Mapfile.hpp"
#include "ssp_shiftprofile.h"

namespace SSP {
  class Global: private Uncopyable {
    bool Greekchr;
    
    std::string samplename;
    std::string oprefix;
    std::string obinprefix;
    bool lackOfRead4FragmentVar;
  
  public:
    SeqStatsGenome genome;
    class LibComp complexity;
    SSPstats sspst;
    std::vector<SeqWigStats>::iterator lchr; // longest chromosome
  
  Global(): Greekchr(false), lackOfRead4FragmentVar(false), complexity() {}
    
    void setOpts(MyOpt::Opts &allopts) {
      genome.setOpts(allopts);
      sspst.setOpts(allopts);
      complexity.setOpts(allopts);
    }
    void setValues(const MyOpt::Variables &values) {
      genome.setValues(values);
      complexity.setValues(values);
      lchr = setlchr(genome);
      samplename = values["output"].as<std::string>();
      oprefix = values["odir"].as<std::string>() + "/" + values["output"].as<std::string>();
    }
    
    std::string getprefix() const { return oprefix; }      
    void outputSSPstats() {
      std::string filename = getprefix() + ".stats.txt";
      std::ofstream out(filename);
      out << "Sample\ttotal read number\tnonredundant read number\t"
	  << "read length\tfragment length\t";
      sspst.printhead(out);
      out << samplename << "\t" << genome.getnread(Strand::BOTH) << "\t" << genome.getnread_nonred(Strand::BOTH) << "\t"
	  << genome.dflen.getlenF3() << "\t" << genome.dflen.getflen() << "\t";
      sspst.print(out);
    }
    
  };
}

#endif /* _SSP_GV_H_ */
