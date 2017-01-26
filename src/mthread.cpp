/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "mthread.h"

namespace MyMthread {
  
  std::vector<chrrange> getVsepchr(const uint64_t genomelen, const std::vector<SeqStats> &chr, const int32_t numthreads)
  {
    std::vector<chrrange> vsep;
    
    uint32_t sepsize = genomelen/numthreads;
    for(uint32_t i=0; i<chr.size(); ++i) {
      uint32_t s = i;
      uint64_t l(0);
      while(l < sepsize && i<chr.size()) {
	l += chr[i].getlen();
	i++;
      }
      i--;
      uint32_t e = i;
      chrrange sep(s,e);
      vsep.push_back(sep);
    }
    return vsep;
  }
  
}
