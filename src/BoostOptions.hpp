/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _BOOSTOPTIONS_HPP_
#define _BOOSTOPTIONS_HPP_

#include <boost/program_options.hpp>

namespace MyOpt {
  using Variables = boost::program_options::variables_map;
  using Opts      = boost::program_options::options_description;
}

#endif /* _BOOSTOPTIONS_HPP_ */
