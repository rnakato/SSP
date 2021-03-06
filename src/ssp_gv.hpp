/* Copyright(c) Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _SSP_GV_HPP_
#define _SSP_GV_HPP_

#include <fstream>
#include "MThread.hpp"
#include "LibraryComplexity.hpp"
#include "Mapfile.hpp"
#include "ShiftProfile.hpp"
#include "FragmentClusterScore.hpp"

namespace SSP {
  class Global: private Uncopyable {

    std::string samplename;
    std::string oprefix;
    std::string obinprefix;

  public:
    SeqStatsGenomeSSP genome;
    LibComp complexity;
    SSPstats sspst;
    FCSstats fcsst;

    Global():
      complexity(),
	      sspst(5*NUM_100K, NUM_1M, 5000, 500, 1500)
    {}

    void setOpts(MyOpt::Opts &allopts) {
      genome.setOpts(allopts);
      sspst.setOpts(allopts);
      fcsst.setOpts(allopts);
      complexity.setOpts(allopts);
    }
    void setValues(const MyOpt::Variables &values) {
      genome.setValues(values);
      complexity.setValues(values);
      sspst.setValues(values);
      fcsst.setValues(values);
      samplename = MyOpt::getVal<std::string>(values, "output");
      oprefix = MyOpt::getVal<std::string>(values, "odir") + "/" + MyOpt::getVal<std::string>(values, "output");
    }

    std::string getprefix() const { return oprefix; }
    void outputSSPstats() {
      std::string filename = getprefix() + ".stats.txt";
      std::ofstream out(filename);
      out << "Sample\ttotal read number\tnonredundant read number\tread length\tfragment length\t";
      sspst.printhead(out);
      out << "\t";
      fcsst.printhead(out);
      out << samplename << "\t"
          << genome.getnread(Strand::BOTH) << "\t"
          << genome.getnread_nonred(Strand::BOTH) << "\t"
          << genome.dflen.getlenF3() << "\t"
          << genome.dflen.getflen() << "\t";
      sspst.print(out);
      out << "\t";
      fcsst.print(out);
    }

  };
}

#endif /* _SSP_GV_HPP_ */
