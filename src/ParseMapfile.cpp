/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <htslib/sam.h>
#include "Mapfile.hpp"
#include "../common/gzstream.h"

namespace {
  void addFragToChr(SeqStatsGenomeSSP &genome, const Fragment &frag)
  {
    genome.dflen.addF3(frag.readlen_F3);
    genome.dflen.addvflen(frag.fraglen);

    bool on(false);
    for (auto &x: genome.chr) {
      if (x.getname() == frag.chr) {
	x.addfrag(frag);
	on = true;
      }
    }
    if (!on) std::cerr << "Warning: " << frag.chr << " is not in genometable." << std::endl;

    return;
  }

  void do_bampe(SeqStatsGenomeSSP &genome, const std::string &inputfile)
  {

    htsFile *fp = hts_open(inputfile.c_str(), "r");  //open bam file
    bam_hdr_t *bamHdr = sam_hdr_read(fp); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
//    printf("%s\n", bamHdr->text);

    uint64_t mappedReads(0), matchedReads(0), matchedProperReads(0), unmatchedReads(0), ProperPair(0);
    uint64_t unmappedReads(0);
    uint64_t duplicatedReads(0);
    uint64_t failqualityReads(0);
    uint64_t forwardReads(0), reverseReads(0);

    while (sam_read1(fp, bamHdr, aln) >= 0) { // 0: SAM, >0: BAM, CRAM
      auto &x = aln->core;
      uint32_t flag = x.flag;
      bool is_mapped = !(flag&4);

      if (is_mapped) {
	++mappedReads;

	std::string chr = bamHdr->target_name[x.tid];
	int32_t position = x.pos;
	int32_t readlen = x.l_qseq;
	int32_t isize  = x.isize;
	bool is_failquality = flag&512;
	bool is_duplicate = flag&1024;
	bool is_paired = flag&1;
	bool is_properpair = flag&2;
	bool is_1stread_of_pair = flag&64;
	bool is_unmatched_pair = flag&8;

	if (is_duplicate) {
	  ++duplicatedReads;
	  continue;
	}
	if (is_failquality) {
	  ++failqualityReads;
	  continue;
	}
	if (!is_paired) continue;

	if (is_unmatched_pair) ++unmatchedReads;
	else {
	  ++matchedReads;
	  if (is_properpair) ++matchedProperReads;
	}
	if (!is_1stread_of_pair) genome.dflen.addF5(readlen);

	if (is_properpair && is_1stread_of_pair) {
	  ++ProperPair;
	  bool strand = bam_is_rev(aln); // 0: forward 1: reverse
	  if (strand) ++reverseReads; else ++forwardReads;

	  Fragment frag;
	  frag.addSAM(chr, readlen, position, isize, strand, genome.isPaired());
	  if (frag.fraglen > genome.getmaxins()) continue;
	  frag.print();
	  addFragToChr(genome, frag);
	}
      }
      else ++unmappedReads;
    }

    bam_destroy1(aln);
    sam_close(fp);

    std::cout << "mapped reads: " << mappedReads
	      << "\tunmapped reads: " << unmappedReads
	      << "\nmatched reads: " << matchedReads
	      << "\tmatched proper reads: " << matchedProperReads
	      << "\tunmatched reads: " << unmatchedReads
	      << "\nproperpair: " << ProperPair
	      << "\t+ pairs: " << forwardReads
	      << "\t- pairs: " << reverseReads
	      << "\nduplicated reads: " << duplicatedReads
	      << "\nFalied quality reads: " << failqualityReads
	      << std::endl;

    return;
  }

  void PrintPairWarning(const bool is_paired)
  {
    static bool pairwarning(false);
    if (is_paired) {
      if (!pairwarning) {
	std::cerr << "Warning: parsing paired-end file as single-end." << std::endl;
	pairwarning = true;
      }
    }
    return;
  }

  void do_bamse(SeqStatsGenomeSSP &genome, const std::string &inputfile)
  {
    htsFile *fp = hts_open(inputfile.c_str(), "r");  //open bam file
    bam_hdr_t *bamHdr = sam_hdr_read(fp); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment
//    printf("%s\n", bamHdr->text);

    uint64_t mappedReads(0);
    uint64_t unmappedReads(0);
    uint64_t duplicatedReads(0);
    uint64_t failqualityReads(0);
    uint64_t forwardReads(0), reverseReads(0);

    while (sam_read1(fp, bamHdr, aln) >= 0) { // 0: SAM, >0: BAM, CRAM
      auto &x = aln->core;
      uint32_t flag = x.flag;
      bool is_mapped = !(flag&4);

      if (is_mapped) {
	++mappedReads;

	std::string chr = bamHdr->target_name[x.tid];
	int32_t position = x.pos;
	int32_t readlen = x.l_qseq;
	//      uint32_t q2  = x.qual ; //mapping quality
	//      uint8_t *q   = bam_get_seq(aln); //quality string
	//      int32_t mpos  = x.mpos;
	int32_t isize  = x.isize;
	bool is_failquality = flag&512;
	bool is_duplicate = flag&1024;
	bool is_paired = flag&1;
	if (is_duplicate) {
	  ++duplicatedReads;
	  continue;
	}
	if (is_failquality) {
	  ++failqualityReads;
	  continue;
	}
	bool strand = bam_is_rev(aln); // 0: forward 1: reverse
	if (strand) ++reverseReads; else ++forwardReads;

	PrintPairWarning(is_paired);

	Fragment frag;
	frag.addSAM(chr, readlen, position, isize, strand, genome.isPaired());
	frag.print();
	addFragToChr(genome, frag);
      }
      else ++unmappedReads;
    }

    bam_destroy1(aln);
    sam_close(fp);

    std::cout << "mapped reads: " << mappedReads
	      << "\t+ reads: " << forwardReads
	      << "\t- reads: " << reverseReads
	      << "\nduplicated reads: " << duplicatedReads
	      << "\nFalied quality reads: " << failqualityReads
	      << "\nunmapped reads: " << unmappedReads
	      << std::endl;

    return;
  }

  void parseSam(const std::string &inputfile, SeqStatsGenomeSSP &genome)
  {
    if ((genome.onFtype() && genome.getftype() == "SAM") || isStr(inputfile, ".sam")) {
      std::cout << "Input format: SAM" << std::endl;
    } else if ((genome.onFtype() && genome.getftype() == "BAM") || isStr(inputfile, ".bam")) {
      std::cout << "Input format: BAM" << std::endl;
    } else if ((genome.onFtype() && genome.getftype() == "CRAM") || isStr(inputfile, ".cram")) {
      std::cout << "Input format: CRAM" << std::endl;
    } else {
      PRINTERR_AND_EXIT("error: invalid input file type.");
    }

    if (genome.isPaired()) do_bampe(genome, inputfile);
    else do_bamse(genome, inputfile);

    return;
  }

  void parseBowtie(const std::string &inputfile, SeqStatsGenomeSSP &genome)
  {
    std::ifstream in(inputfile);
    if (!in) PRINTERR_AND_EXIT("Could not open " << inputfile << ".");
    std::cout << "Input format: BOWTIE" << std::endl;

    std::string chr_F3(""), chr_F5(""), nametemp("");
    int32_t F5(0);
    Fragment fragpair;

    std::string lineStr;
    while (!in.eof()) {
      getline(in, lineStr);
      if (lineStr.empty()) continue;

      std::vector<std::string> v;
      ParseLine(v, lineStr, '\t');
      //      ParseLine(v, lineStr, "\t");
      //      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

      if (genome.isPaired()) {
	std::vector<std::string> read;
	ParseLine(read, v[0], '/');
	//	ParseLine(read, v[0], "/");
	//boost::split(read, v[0], boost::algorithm::is_any_of("/"));
	if (nametemp != "" && nametemp != read[0]) PRINTERR_AND_EXIT("Invalid read pair." << nametemp <<"-" << read[0]);
	if (read[1] == "1") {  // F3 read
	  chr_F3 = rmchr(v[2]);
	  fragpair.readlen_F3 = v[4].length();
	  if (v[1] == "+") {
	    fragpair.strand = Strand::FWD;
	    fragpair.F3 = stoi(v[3]);
	  } else {
	    fragpair.strand = Strand::REV;
	    fragpair.F3 = stoi(v[3]) + fragpair.readlen_F3;
	  }
	} else {
	  chr_F5 = rmchr(v[2]);
	  if (v[1] == "+") F5 = stoi(v[3]);
	  else            F5 = stoi(v[3]) + v[4].length();
	  genome.dflen.addF5(v[4].length());
	}
	if (chr_F3 != "" && chr_F5 != ""){
	  if (chr_F3 == chr_F5) {
	    fragpair.chr = chr_F3;
	    fragpair.fraglen = abs(F5 - fragpair.F3);
	    if (fragpair.fraglen <= genome.getmaxins()) addFragToChr(genome, fragpair);
	    fragpair.print();
	  }
	  chr_F3 = "";
	  chr_F5 = "";
	  nametemp = "";
	}
      } else {
	if (isStr(v[0], "/2")) PRINTERR_AND_EXIT("Warning: parsing paired-end file as single-end");
	Fragment frag;
	frag.chr = rmchr(v[2]);
	frag.readlen_F3 = v[4].length();
	if (v[1] == "+") {
	  frag.strand = Strand::FWD;
	  frag.F3 = stoi(v[3]);
	} else {
	  frag.strand = Strand::REV;
	  frag.F3 = stoi(v[3]) + frag.readlen_F3;
	}
	frag.print();
	addFragToChr(genome, frag);
      }

    }
    return;
  }

  template <class T>
  void funcTagAlign(SeqStatsGenomeSSP &genome, T &in)
  {
    std::cout << "Input format: TAGALIGN" << std::endl;
    std::string lineStr;
    while (!in.eof()) {
      getline(in, lineStr);
      if (lineStr.empty()) continue;

      std::vector<std::string> v;
      ParseLine(v, lineStr, '\t');
      //      ParseLine(v, lineStr, "\t");
      //      boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));
      if (v.size() < 6) PRINTERR_AND_EXIT("Use tagAlign (BED3+3) file");

      if (genome.isPaired()) PRINTERR_AND_EXIT("tagAlign format does not support paired-end file.\n");
      else {
	int32_t start(stoi(v[1]));
	int32_t end(stoi(v[2]));
	Fragment frag;
	frag.chr = rmchr(v[0]);
	frag.readlen_F3 = abs(end - start);
	if (v[5] == "+") {
	  frag.strand = Strand::FWD;
	  frag.F3 = start;
	} else {
	  frag.strand = Strand::REV;
	  frag.F3 = start + frag.readlen_F3;
	}
	frag.print();
	addFragToChr(genome, frag);
      }
    }
    return;
  }

  void parseTagAlign(const std::string &inputfile, SeqStatsGenomeSSP &genome)
  {
    if (isStr(inputfile, ".gz")) {

      igzstream in(inputfile.c_str());
      funcTagAlign(genome, in);
    } else {
      std::ifstream in(inputfile);
      if (!in) PRINTERR_AND_EXIT("Could not open " << inputfile << ".");
      funcTagAlign(genome, in);
    }
    return;
  }

  /*  void parseTagAlign(const std::string &inputfile, SeqStatsGenomeSSP &genome)
  {
    if (isStr(inputfile, ".gz")) {
      std::string command = "zcat " + inputfile;
      FILE *fp = popen(command.c_str(), "r");
      __gnu_cxx::stdio_filebuf<char> *p_fb = new __gnu_cxx::stdio_filebuf<char>(fp, std::ios_base::in);
      std::istream in(static_cast<std::streambuf *>(p_fb));
      funcTagAlign(genome, in);
    } else {
      std::ifstream in(inputfile);
      if (!in) PRINTERR_AND_EXIT("Could not open " << inputfile << ".");
      funcTagAlign(genome, in);
    }
    return;
    }*/

  /*int32_t check_sv(int32_t sv)
{
  // for paired-end
   LOG("   the read is paired in sequencing: %d\n",sv&1);
  LOG("   the read is mapped in a proper pair: %d\n",sv&2);
  LOG("   the query sequence itself is unmapped: %d\n",sv&4);
  LOG("   the mate is unmapped: %d\n",sv&8);
  LOG("   strand of the query (1 for reverse): %d\n",sv&16);
  LOG("   strand of the mate: %d\n",sv&32);
  LOG("   the read is the first read(F3) in a pair: %d\n",sv&64);
  LOG("   the read is the second read(F5) in a pair: %d\n",sv&128);
  LOG("   the alignment is not primary: %d\n",sv&256);
  LOG("   the read fails platform/vendor quality checks: %d\n",sv&512);
  LOG("   the read is either a PCR or an optical duplicate: %d\n",sv&1024);

  LOG("   template having multiple segments in sequencing: %d\n",sv&1);
  LOG("   each segment properly aligned according to the aligner: %d\n",sv&2);
  LOG("   segment unmapped: %d\n",sv&4);
  LOG("   next segment in the template unmapped: %d\n",sv&8);
  LOG("   SEQ being reverse complemented: %d\n",sv&16);
  LOG("   SEQ of the next segment in the template being reversed: %d\n",sv&32);
  LOG("   the first segment in the template: %d\n",sv&64);
  LOG("   the last segment in the template: %d\n",sv&128);
  LOG("   secondary alignment: %d\n",sv&256);
  LOG("   not passing quality controls: %d\n",sv&512);
  LOG("   PCR or optical duplicate: %d\n",sv&1024);
  LOG("   supplementary alignment: %d\n",sv&2048);


  // unmapped reads
  if (sv&4) goto err;
  // low quality reads
  if (sv&512 || sv&1024) goto err;
  //  if (p->rtype==READTYPE_PAIR){
    // unproper pair
    if (!(sv&2)) goto err;
    // unmatched pairs and interchromosomal pairs
    if (sv&8) goto err;
    // read pair mapped in same strand (for paired-end)
    if ((sv&16 && sv&32) || (!(sv&16) && !(sv&32))) goto err;
    // }

 err:
  return 0;
}  */

}

void SeqStatsGenomeSSP::read_mapfile()
{
  std::vector<std::string> v;
  ParseLine(v, getInputfile(), ',');

  for (auto &inputfile: v) {
    isFile(inputfile);
    std::cout << boost::format("Parsing %1%...\n") % inputfile;
    if (onFtype()) {
      if (getftype() == "SAM"  || getftype() == "BAM" || getftype() == "CRAM") parseSam(inputfile, *this);
      else if (getftype() == "BOWTIE")   parseBowtie(inputfile, *this);
      else if (getftype() == "TAGALIGN") parseTagAlign(inputfile, *this);
    } else {
      if (isStr(inputfile, ".sam") || isStr(inputfile, ".bam") || isStr(inputfile, ".cram")) parseSam(inputfile, *this);
      else if (isStr(inputfile, ".bowtie"))     parseBowtie(inputfile, *this);
      else if (isStr(inputfile, ".tagalign"))   parseTagAlign(inputfile, *this);
    }
  }

  if (!getnread(Strand::BOTH)) {
    if (isPaired()) PRINTERR_AND_EXIT("No read in input file. Is this a single-end reads file?");
    else PRINTERR_AND_EXIT("No read in input file.");
  }

  return;
}
