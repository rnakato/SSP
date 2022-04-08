# Changelog

## 1.2.3 (2022-04-09)
- Bug fix: fixed the error when the F3 length (read length) is too short (e.g., 1~3)

## 1.2.2 (2020-10-06)
- Bug fix: switch boost::bind to std::bind to avoid complilation error depend on the version of compiler

## 1.2.1 (2020-05-25)
- Adopt CMake to compile SSP
- Add script change_chrname_[to|from]_Greek.pl
- Remove common/statistics.*

## 1.2.0 (2020-05-19)
- Adopt htslib API for parsing SAM/BAM/CRAM format
- Now available to use on Mac
- Remove BedFormat.hpp that is used in DROMPAplus
- Remove unused classes in SSP that is for DROMPAplus

## 1.1.6 (2020-05-17)
- Change C++ compiler from g++ to clang++

## 1.1.5 (2020-03-12)
- Modify Boost::options to get more intelligible error messages

## 1.1.4 (2020-02-25)
- Remove unnecesary source files from common/ directory

## 1.1.3 (2019-04-15)
- Accept CRAM format for the input file
- script/makegenometable.pl: accept scaffold in Ensembl genome

## 1.1.2 (2018-07-22)
- Chenge the prefix of ReadLengthDist.csv and FragmentLengthDist.csv to *.tsv
- Add --allchr options (for DROMPAplus)

## 1.1.1 (2018-04-25)
- Improve to read BAM file faster

## 1.1.0 (2018-02-10)
- Bug fix in --ssp_cc option
- Add genometable files in data directory

## 1.0.2 (2017-09-04)
- Add LICENSE
- Add ChangeLog
- Remove alglib library from source

## 1.0.1 (2017-07-19)
- Bug fix in --pair option

## 1.0.0
- First commit
