/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#include "ssp_estFlen.h"
#include "ssp_shiftprofile.h"
#include <time.h>

void estimateFragLength(const MyOpt::Variables &values, Mapfile &p)
{
  //  if(!values.count("pair") && !values.count("nomodel")) {
    clock_t t1 = clock();
    strShiftProfile(values, p, "exjaccard");
    clock_t t2 = clock();
    std::cout << "Jaccard Vec: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";

    strShiftProfile(values, p, "jaccard"); 
    clock_t t3 = clock();
    std::cout << "Jaccard Bit: " << static_cast<double>(t3 - t2) / CLOCKS_PER_SEC << "sec.\n";

    strShiftProfile(values, p, "fvp"); 
    clock_t t4 = clock();
    std::cout << "Fragment variability: " << static_cast<double>(t4 - t3) / CLOCKS_PER_SEC << "sec.\n";

    strShiftProfile(values, p, "hdp");
    clock_t t5 = clock();
    std::cout << "Hamming: " << static_cast<double>(t5 - t4) / CLOCKS_PER_SEC << "sec.\n";

    strShiftProfile(values, p, "ccp");
    clock_t t6 = clock();
    std::cout << "ccp: " << static_cast<double>(t6 - t5) / CLOCKS_PER_SEC << "sec.\n";
    //  }
}
