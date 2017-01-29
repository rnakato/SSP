/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#ifndef _PARSEBED_HPP_
#define _PARSEBED_HPP_

template <class T>
std::vector<T> parseBed(const std::string &fileName)
{
  std::vector<T> vbed;
  std::ifstream in(fileName);
  if(!in) PRINTERR("BED file does not exist.");

  std::string lineStr;
  std::vector<std::string> v;
  while (!in.eof()) {
    getline(in, lineStr);
    
    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
    if(v[1] == "start") continue;
    T bed(v);
    vbed.push_back(bed);
  }

  return vbed;
}

template <class T>
void printBed(const std::vector<T> &vbed)
{
  for (auto &x: vbed) {
    x.print();
    std::cout << std::endl;
  }
  std::cout << "bed num: " << vbed.size() << std::endl;
  return;
}


#endif  // _PARSEBED_HPP_
