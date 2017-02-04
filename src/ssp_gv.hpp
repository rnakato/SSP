/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SSP_GV_HPP_
#define _SSP_GV_HPP_

#include <fstream>
#include "MThread.hpp"
#include "LibraryComplexity.hpp"
#include "Mapfile.hpp"
#include "ssp_shiftprofile.hpp"

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
  
  Global(): Greekchr(false), lackOfRead4FragmentVar(false), complexity() {}
    
    void setOpts(MyOpt::Opts &allopts) {
      genome.setOpts(allopts);
      sspst.setOpts(allopts);
      complexity.setOpts(allopts);
    }
    void setValues(const MyOpt::Variables &values) {
      genome.setValues(values);
      complexity.setValues(values);
      sspst.setValues(values);
      samplename = MyOpt::getVal<std::string>(values, "output");
      oprefix = MyOpt::getVal<std::string>(values, "odir") + "/" + MyOpt::getVal<std::string>(values, "output");
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

#endif /* _SSP_GV_HPP_ */
