/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "util.hpp"
#include <algorithm>
#include <boost/filesystem.hpp>

void ParseLine(std::vector<std::string> &v, const std::string &str, char delim)
{
  size_t current(0), found;
  while((found = str.find_first_of(delim, current)) != std::string::npos) {
    v.emplace_back(std::string(str, current, found - current));
    current = found + 1;
  }
  v.emplace_back(std::string(str, current, str.size() - current));
  return;
}
  
void printList()
{
  std::cout << std::endl;
}

std::string rmchr(const std::string &chr)
{
  std::string s;
  if(!chr.find("chr")) s = chr.substr(3);
  else s = chr;
  return s;
}

bool checkFile(const std::string &str)
{
  boost::filesystem::path const file(str);
  return boost::filesystem::exists(file);
}

void isFile(const std::string &str)
{
  boost::filesystem::path const file(str);
  if(!boost::filesystem::exists(file)) {
    std::cerr << "Error: " << str << " does not exist." << std::endl;
    exit(1);
  }
}

bool isStr(std::string str, std::string query)
{
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::transform(query.begin(), query.end(), query.begin(), ::tolower);
  if(str.find(query) != std::string::npos) return true;
  else return false;
}
