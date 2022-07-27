/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include "util.hpp"
#include <algorithm>
#include <boost/filesystem.hpp>

void ParseLine_NoDelimCheck(std::vector<std::string> &v, const std::string &str, char delim)
{
  size_t current(0), found(0);

  while((found = str.find_first_of(delim, current)) != std::string::npos) {
    v.emplace_back(std::string(str, current, found - current));
    current = found + 1;
  }
  v.emplace_back(std::string(str, current, str.size() - current));

  return;
}

int32_t ParseLine(std::vector<std::string> &v, const std::string &str, char delim)
{
  size_t current(0), found(0);

  if(str.find(delim) == std::string::npos) {    // no delim in str
    return 1;
  }

  while((found = str.find_first_of(delim, current)) != std::string::npos) {
    v.emplace_back(std::string(str, current, found - current));
    current = found + 1;
  }
  v.emplace_back(std::string(str, current, str.size() - current));

  return 0;
}

void printList()
{
  std::cout << std::endl;
}

std::string rmchr(const std::string &chr)
{
  std::string s;
  if(!chr.find("chr") || !chr.find("Chr")) s = chr.substr(3);
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

int32_t compare_chr(const std::string &chr1, const std::string &chr2) {
  int32_t c1, c2;
  try {
    c1 = stoi(chr1);
    c2 = stoi(chr2);
    return c1 - c2;
  } catch (const std::invalid_argument& e) {
    try {
      c1 = stoi(chr1);
      return -1; // c2の方が大きい
    } catch (const std::invalid_argument& e) {
      try {
	c2 = stoi(chr2);
	return 1; // c1の方が大きい
      } catch (const std::invalid_argument& e) {
	if (chr1 == chr2) return 0;
	else if (chr1 == "X") return -1;
	else if (chr1 == "Y") {
	  if (chr2 == "X") return 1; else return -1;
	}
	else return 1;
      }
    }
  }
  return 0;
}
