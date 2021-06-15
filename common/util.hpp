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

void ParseLine(std::vector<std::string> &v, const std::string &str, char delim);
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

/*template <class T>
void GaussianSmoothing(std::vector<T> &v, const int32_t nsmooth)
{
  std::vector<double> w(nsmooth+1,0);
  double var(1);

  double sum(0);
  for (int32_t i=0; i<=nsmooth; ++i) {
    w[i] = exp(static_cast<double>(-i*i)/2*var*var);
    sum += w[i];
  }
  double r(1/(sum*2 - w[0]));

  std::vector<double> m(nsmooth+1,0);
  for (int32_t i=0; i<=nsmooth; ++i) m[i] = v[nsmooth-i];

  //  std::cout << "GS Weight: ";
  //for (int32_t i=0; i<=nsmooth; ++i)  std::cout << w[i] << "\t";
  //std::cout  << "r: " << r << std::endl;

  for (size_t i=nsmooth; i<v.size()-nsmooth; ++i) {
    m[0] = v[i];

    //    std::cout << "before: ";
    //    for (int32_t j=0; j<nsmooth; ++j) std::cout << v[i-nsmooth+j] << "\t";
    // for (int32_t j=nsmooth; j>=0; --j) std::cout << m[j] << "\t";
    //for (int32_t j=1; j<=nsmooth; ++j) std::cout << v[i+j] << "\t";
    //std::cout << std::endl;
    double val(w[0]*m[0]);
    for (int32_t j=1; j<=nsmooth; ++j) val += w[j] * (m[j] + v[i+j]);
    v[i] = val*r;
    //std::cout << v[i] << std::endl;

    for (int32_t i=nsmooth; i>0; --i) m[i] = m[i-1];
  }
  return;
}*/

template <class T>
int32_t findIndex(std::vector<T> array, T value)
{
    auto iter = std::find(array.begin(), array.end(), value);
    size_t index = std::distance(array.begin(), iter);
    if(index == array.size()) index = -1;

    return index;
}

#endif /* _UTIL_HPP_ */
