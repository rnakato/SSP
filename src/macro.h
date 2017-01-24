/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _MACRO_H_
#define _MACRO_H_

#define VALUE2WIGARRAY(v) ((v) * 1000.0)
#define WIGARRAY2VALUE(v) ((v) / 1000.0)
#define BPRINT std::cout << boost::format
#define PRINTERR(...) do{ std::cerr << "Error: " << __VA_ARGS__ << std::endl; std::exit(1); }while(0)

enum {NUM_1K=1000,
      NUM_100K=100000,
      NUM_1M=1000000,
      NUM_10M=10000000,
      NUM_100M=100000000};

enum PWfile_Type {
  TYPE_BINARY,
  TYPE_COMPRESSWIG,
  TYPE_UNCOMPRESSWIG,
  TYPE_BEDGRAPH,
  TYPE_BIGWIG,
  PWFILETYPENUM
};

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

#endif /* _MACRO_H_ */
