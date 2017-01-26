/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _MTHREAD_H_
#define _MTHREAD_H_

#include <cstdint>
#include <vector>
#include "mapfileclass.h"

namespace MyMthread {
  class chrrange {
  public:
    uint32_t s;
    uint32_t e;
  chrrange(uint32_t start, uint32_t end): s(start), e(end) {}
  };
  
  std::vector<chrrange> getVsepchr(const uint64_t genomelen, const std::vector<SeqStats> &chr, const int32_t numthreads);
}

#endif /* _MTHREAD_H_ */
