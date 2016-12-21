/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of SSP sources.
 */
#ifndef _SSP_SHIFTPROFILE_H_
#define _SSP_SHIFTPROFILE_H_

#include "pw_gv.h"

void strShiftProfile(const MyOpt::Variables &values, Mapfile &p, std::string);
void makeFCSProfile(const MyOpt::Variables &values, Mapfile &p, const std::string &typestr);

#endif /* _SSP_SHIFTPROFILE_H_ */
