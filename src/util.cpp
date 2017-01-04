/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "util.h"
#include <algorithm>
#include <boost/filesystem.hpp>

std::string rmchr(const std::string &chr)
{
  std::string s;
  if(!chr.find("chr")) s = chr.substr(3);
  else s = chr;
  return s;
}

void isFile(const std::string &str)
{
  boost::filesystem::path const file(str);
  if(!boost::filesystem::exists(file)) {
    std::cerr << "Error: " << str << " does not exist." << std::endl;
    exit(1);
  }
}

int isStr(std::string str, std::string query)
{
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::transform(query.begin(), query.end(), query.begin(), ::tolower);
  if(str.find(query) != std::string::npos) return 1;
  else return 0;
}
