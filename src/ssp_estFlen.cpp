/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#include "ssp_estFlen.h"
#include "ssp_shiftprofile.h"
#include <time.h>

void estimateFragLength(const MyOpt::Variables &values, Mapfile &p)
{
  if(values.count("pair") || values.count("nomodel")) return;

  clock_t t1,t2;
  t1 = clock();
  strShiftProfile(values, p, "jaccard"); 
  t2 = clock();
  std::cout << "Jaccard Bit: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";

  strShiftProfile(values, p, "fvp"); 
  clock_t t3 = clock();
  std::cout << "Fragment variability: " << static_cast<double>(t3 - t2) / CLOCKS_PER_SEC << "sec.\n";

  if(values.count("ssp_exjac")) {
    t1 = clock();
    strShiftProfile(values, p, "exjaccard");
    t2 = clock();
    std::cout << "Jaccard Vec: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }

  if(values.count("ssp_hd")) {
    t1 = clock();
    strShiftProfile(values, p, "hdp");
    t2 = clock();
    std::cout << "Hamming: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }
    
  if(values.count("ssp_cc")) {
    t1 = clock();
    strShiftProfile(values, p, "ccp");
    t2 = clock();    
    std::cout << "ccp: " << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
  }
  
  return;
}
