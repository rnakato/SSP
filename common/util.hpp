/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _UTIL_HPP_
#define _UTIL_HPP_
#include <iostream>
#include <string>
#include <vector>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "inline.hpp"

std::string rmchr(const std::string &chr);
void isFile(const std::string &);
bool checkFile(const std::string &str);
bool isStr(std::string, std::string);
int32_t compare_chr(const std::string &, const std::string &);

void ParseLine_NoDelimCheck(std::vector<std::string> &v, const std::string &str, char delim);
int32_t ParseLine(std::vector<std::string> &v, const std::string &str, char delim);
void printList();

template <class Thead, class... Tbody>
void printList(Thead head, Tbody... body)
{
  std::cout << head;
  if(sizeof...(body) > 0) std::cout << '\t';
  printList(body...);
}

template <class T>
int32_t getmaxi(std::vector<T> v)
{
  T max(0);
  int32_t maxi(0);
  int32_t size = v.size();
  for(int32_t i=0; i<size; ++i) {
    if(max < v[i]) {
      max = v[i];
      maxi = i;
    }
  }
  return maxi;
};

template <class T>
int32_t findIndex(std::vector<T> array, T value)
{
    auto iter = std::find(array.begin(), array.end(), value);
    size_t index = std::distance(array.begin(), iter);
    if(index == array.size()) index = -1;

    return index;
}

#endif /* _UTIL_HPP_ */
