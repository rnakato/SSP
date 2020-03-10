/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _INLINE_H_
#define _INLINE_H_

#include <iostream>
#include <sstream>
#include <boost/format.hpp>
#include "statistics.hpp"

#define PRINTERR_AND_EXIT(...) do{ std::cerr << "Error: " << __VA_ARGS__ << std::endl; std::exit(1); } while(0)
#ifdef DEBUG
#define DEBUGprint(...) do{ std::cout << __VA_ARGS__ << std::endl; } while(0)
#define DEBUGprint_FUNCStart() do{ std::cout << "FUNCstart: " << __PRETTY_FUNCTION__ << std::endl; } while(0)
#define DEBUGprint_FUNCend() do{ std::cout << "FUNCend: " << __PRETTY_FUNCTION__  << std::endl; } while(0)
#else
#define DEBUGprint(...) do{} while(0)
#define DEBUGprint_FUNCStart(...) do{} while(0)
#define DEBUGprint_FUNCend(...) do{} while(0)
#endif

enum {NUM_1K=1000,
      NUM_100K=100000,
      NUM_1M=1000000,
      NUM_10M=10000000,
      NUM_100M=100000000};

class Uncopyable {
 protected:
  Uncopyable(){}
  ~Uncopyable(){}
 private:
  Uncopyable(const Uncopyable &);
  Uncopyable& operator=(const Uncopyable &);
};

template <class T, class S>
inline bool my_range(const T i, const S min, const S max)
{
  bool b;
  if((i >=min) && (i <=max)) b=1; else b=0;
  return b;
}

template <class T>
inline bool my_overlap(const T s1, const T e1, const T s2, const T e2)
{
  return (e1 >= s2) && (e2 >= s1);
}

template <class T, class S>
inline double getratio(const T x, const S y)
{
  //  if(!y) std::cerr << "Warning: denominator=0." << std::endl;
  return y ? x/static_cast<double>(y): 0;
}

template <class T, class S>
inline double getpercent(const T x, const S y)
{
  return getratio(x,y)*100;
}

template <class T, class S>
  inline void printNumandPer(std::ofstream &out, T a, S b)
{
  out << boost::format("%1% (%2$.1f%%)\t") % a % getpercent(a,b);
};

inline void PrintTime(const clock_t t1, const clock_t t2, const std::string &str)
{
#ifdef CLOCK
  std::cout << str << ": "
	    << static_cast<double>(t2 - t1) / CLOCKS_PER_SEC
	    << "sec.\n";
#else
  if(str==""){}
  if(t1==0){}
  if(t2==0){}
#endif
}

template<class A, size_t N, class T>
void Fill(A (&array)[N], const T &val) {
  std::fill( (T*)array, (T*)(array+N), val );
}

#endif /* _INLINE_H_ */
