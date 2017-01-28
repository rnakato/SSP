/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _BOOSTOPTIONS_HPP_
#define _BOOSTOPTIONS_HPP_

#include <iostream>
#include <boost/program_options.hpp>

namespace MyOpt {
  using Variables = boost::program_options::variables_map;
  using Opts      = boost::program_options::options_description;
  using Err       = boost::program_options::validation_error;

  template<typename T>
  void range(const T& v, const T& min, const T& max, const char * optname)
  {
    if (v<min || v>max)
      throw Err	(Err::invalid_option_value, optname, std::to_string(v));
  }
  
  template<typename T>
  void over(const T& v, const T& min, const char * optname)
  {
    if (v < min) 
      throw Err	(Err::invalid_option_value, optname, std::to_string(v));
  }

  void setOptIO(Opts &, const std::string &);
  void setOptPair(Opts &);
  void setOptOther(Opts &);
}

#endif /* _BOOSTOPTIONS_HPP_ */
