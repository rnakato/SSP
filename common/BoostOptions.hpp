/* Copyright(c) Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
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

  template <typename T>
  void range(const T &v, const T &min, const T &max, const char *optname)
  {
    if (v < min || v > max)
      throw Err(Err::invalid_option_value, optname, std::to_string(v));
  }

  template <typename T>
  void over(const T& v, const T& min, const char *optname)
  {
    if (v < min)
      throw Err(Err::invalid_option_value, optname, std::to_string(v));
  }

  template <typename T>
  T getVal(const MyOpt::Variables &values, const std::string &optname)
  {
    try {
#ifdef DEBUG
      std::cout << optname << "," << values.count(optname) << std::endl;
#endif
      return values[optname].as<T>();
    }
    catch (const boost::bad_any_cast& e) {
      std::cerr << "option --" << optname << ": "<< e.what() << std::endl;
      auto str = values[optname].as<std::string>();
      if (str.find("-") == 0) {
	std::cerr << "Error: specify a value to --" << optname << "." << std::endl;
      }
      exit(0);
    }
  }

  template <>
  inline std::string getVal<std::string>(const MyOpt::Variables &values, const std::string &optname)
  {
    try {
#ifdef DEBUG
      std::cout << optname << "," << values.count(optname) << std::endl;
#endif
      auto str = values[optname].as<std::string>();
      if (str.find("-") == 0) {
	std::cerr << "Error: specify a value to --" << optname << "." << std::endl;
	exit(0);
      }
      return values[optname].as<std::string>();
    }
    catch (const boost::bad_any_cast& e) {
      std::cerr << "option --" << optname << ": "<< e.what() << std::endl;
      exit(0);
    }
  }

  template <>
  inline std::vector<std::string> getVal<std::vector<std::string>>(const MyOpt::Variables &values, const std::string &optname)
  {
#ifdef DEBUG
    std::cout << optname << ",tokushuka," << values.count(optname) << std::endl;
#endif
    try {
      auto v = values[optname].as<std::vector<std::string>>();

      int32_t n(1);
      for (auto &x: v) {
	if (x.find("-") == 0) {
	  std::cerr << "Error: specify a value to " << std::to_string(n) << "th --" << optname << "." << std::endl;
	  exit(0);
	}
	++n;
      }
      return v;
    }
    catch (const boost::bad_any_cast& e) {
      std::cerr << "option --" << optname << ": "<< e.what() << std::endl;
      exit(0);
    }
  }

  template <class T>
  void printOpt(const boost::program_options::variables_map &values,
		const std::string &opt, const std::string &str)
  {
    if (values.count(opt)) std::cout << str << ": " << values[opt].as<T>() << std::endl;
    return;
  }

  template <class T>
  void printVOpt(const boost::program_options::variables_map &values,
		 const std::string &opt, const std::string &str)
  {
    if (values.count(opt)) {
      auto v = values[opt].as<std::vector<T>>();
      for(uint32_t i=0; i<v.size(); ++i) {
	std::cout << str << " " << (i+1) << ": " << v[i] << std::endl;
      }
    }
    return;
  }

  void setOptIO(Opts &, const std::string &);
  void setOptPair(Opts &);
  void setOptOther(Opts &);
  void dumpIO(const Variables &);
  void dumpGenomeTable(const Variables &);
  void dumpPair(const Variables &);
  void dumpLibComp(const Variables &);
  void dumpFragmentLengthDist(const Variables &);
  void dumpOther(const Variables &);
}

#endif /* _BOOSTOPTIONS_HPP_ */
