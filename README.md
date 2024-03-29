# SSP README

## 0. Changelog
See [Changelog](https://github.com/rnakato/SSP/blob/master/ChangeLog.md)

# 1. Overview
SSP (strand-shift profile) is a tool for quality assessment of ChIP-seq data without peak calling.


SSP provides metrics to:
- quantify the S/N for both point- and broad-source factors (NSC),
- estimate peak reliability based on the mapped-read distribution throughout a genome (Bu),
- and estimate peak intensity and peak mode (point- or broad-source, FCS).

The outputs of SSP are displayed in PDF format and also written to text files.

# 2. Install

## 2.1 Dependencies
SSP is written in C++11 and requires the following programs and libraries.
The version numbers listed have been tested successfully.
* [Boost C++ library (1.53.0, 1.58.0)](http://www.boost.org/)
* [GNU Scientific Library (1.15, 2.1)](http://www.gnu.org/software/gsl/)
* [zlib (1.2.7, 1.2.8)](http://www.zlib.net/)
* [CMake (>2.8)](https://cmake.org/)
* [HTSlib (1.10.2)](https://github.com/samtools/htslib) (for SAM/BAM/CRAM formatted input)

## 2.2 Building from source

### 2.2.1 Install required libraries

On Ubuntu:

    sudo apt install git build-essential libboost-all-dev libcurl4-gnutls-dev \
                     liblzma-dev libz-dev libbz2-dev cmake

On CentOS and Red Hat:

    sudo yum -y install git gcc-c++ clang boost-devel zlib-devel bzip2-devel cmake

On Mac:

    brew install curl xz zlib boost cmake

### 2.2.2 Install SSP

    git clone https://github.com/rnakato/SSP.git
    cd SSP
    make

### 2.2.3 Add the PATH environment variable
For example, if you downloaded SSP into the $HOME/my_chipseq_exp directory, type:

    export PATH = $PATH:$HOME/my_chipseq_exp/SSP/bin

## 2.3 Docker image

SSP and DROMPA are also probatively available on Docker Hub.

To obtain a docker image for SSP and DROMPA, type:

    docker pull rnakato/ssp_drompa
    docker run -it --rm rnakato/ssp_drompa ssp

For Singularity:

    singularity build ssp_drompa.img docker://rnakato/ssp_drompa
    singularity exec ssp_drompa.img ssp

# 3. Usage
## 3.1. Options
    Usage: ssp [option] -i <inputfile> -o <output> --gt <genome_table>

    Options:

    Input/Output:
      -i [ --input ] arg           Mapping file. Multiple files are allowed (separated by ',')
      -o [ --output ] arg          Prefix of output files
      --odir arg (=sspout)         output directory name
      -f [ --ftype ] arg           {SAM|BAM|CRAM|BOWTIE|TAGALIGN}: format of input file
                                    TAGALIGN can be gzip'ed (extension: tagAlign.gz)

    For paired-end:
      --pair                       add when the input file is paired-end
      --maxins arg (=500)          maximum fragment length

    Genome:
      --gt arg                     Genome table (tab-delimited file describing the name and length of
                                   each chromosome)
      --mptable arg                Genome table of mappable regions
      --include_allchr             Include all chromosomes for calculation (default: autosomes only,
                                   i.e., 'chrN', where N is a numeric number)
    Fragment:
      --nomodel                    omit fraglent length estimation (default: estimated by strand-shift profile)
      --flen arg (=150)            predefined fragment length (with --nomodel option)

    Strand shift profile:
      --num4ssp arg (=10000000)    Read number for calculating backgroud uniformity (per 100 Mbp)
      --ng_from arg (=500000)      start shift of background
      --ng_to arg (=1000000)       end shift of background
      --ng_step arg (=5000)        step shift on of background
      --ssp_cc                     make ssp based on cross correlation
      --ssp_hd                     make ssp based on hamming distance
      --ssp_exjac                  make ssp based on extended Jaccard index
      --eachchr                    make chromosome-sparated ssp files

    Fragment cluster score:
      --ng_from_fcs arg (=100000)  fcs start of background
      --ng_to_fcs arg (=1000000)   fcs end of background
      --ng_step_fcs arg (=100000)  fcs step on of background

    Library complexity:
      --thre_pb arg (=0)           PCRbias threshold (default: more than max(1 read, 10 times greater
                                   than genome average))
      --ncmp arg (=10000000)       read number for calculating library complexity
      --nofilter                   do not filter PCR bias

    Others:
      -p [ --threads ] arg (=1)    number of threads to launch
      -v [ --version ]             print version
      -h [ --help ]                show help message



## 3.2. Tutorial
The simplest command is:

    ssp -i ChIP.sam -o ChIP --gt genometable.txt
then the output files (prefix: "ChIP") are generated in the directory "sspout (default)".
The format of input file is automatically detected by postfix(.sam/.bam/.cram/.bowtie/.tagalign(.gz)). If the detection does not work well, supply -f option (e.g., "-f BAM").

The genome table file (genometable.txt) is a tab-delimited file describing the name and length of each chromosome (see 4.1.)
The chromosome names in the map file and the genome table file must be same.

To supply the mappable genome table and use multiple CPUs:

     ssp -i ChIP.bam -o ChIP --gt genometable.txt --mptable mptable.txt -p 4
"-p 4" specifies the number of CPUs used. The mappable genome table file is necessary for accurate estimation of background uniformity.


SSP allows multiple input files (separated by ",")

     ssp -i ChIP1.bam,ChIP2.bam,ChIP3.bam -o ChIP --gt genometable.txt


Note that the chromosome length should be enough longer than the background length specified. For small genomes (e.g., yeast), the background length should be shorten:

     ssp -i ChIP1.bam -o ChIP --gt genometable.txt --ng_from 10000 --ng_to 50000 --ng_step 500

In this parameter set, the background region is the average ranging from 10k to 50k at steps of 500 bp.

By default, FCS is calcutated for 10M nonredundant reads. If the number of nonredundant reads in the input data are smaller than 10M, specify smaller number for fair comparison among samples as follows:

     ssp -i ChIP1.bam -o ChIP --gt genometable.txt --num4ssp 5000000

When specifying smaller read number for --num4ssp, FCS score becomes smaller, but the magnitude relation among samples is consistent.

By default, SSP uses only autosomes, i.e., ‘chrN’, where N is a numeric number, in the genome_table file. However, in some cases the chromosome names do not start with “chr”. In such a case, add the `--include_allchr` option to include all chromosomes in the genome_table file.

     ssp -i ChIP1.bam -o ChIP --gt genometable.txt --include_allchr

## 3.3. Output files
* ChIP.stats.txt: Stats of the sample (read number, read length, estimated fragment length, NSC, RLSC, RSC, background uniformity, FCS)
* ChIP.jaccard.csv: Jaccard score for each strand shift d
* ChIP.jaccard.pdf: Strand-shift profiles (-500 < d < 1500 and 0 < d < 1M)
* ChIP.jaccard.R: R script to make ChIP.jaccard.pdf
* ChIP.jaccard.R.log: Log file of ChIP.jaccard.R
* ChIP.fcs.csv: FCS for each strand shift d
* ChIP.pnf.csv: PNF (the proportion of neighboring fragments) for each s
* ChIP.FCS.pdf: Profiles of PNF, cPNF (the cumulative proportion of neighboring fragments) and FCS
* ChIP.FCS.R, ChIP.FCS.R.log: R script and log file to make ChIP.FCS.pdf

# 4. Annotation files
## 4.1. Genome table
The genome table file is a tab-delimited file describing the name and length of each chromosome.
To make it, use makegenometable.pl in scripts directory as follows:

     scripts/makegenometable.pl genome.fa > genometable.txt

## 4.2. Mappability table
The mappability table file is a tab-delimited file describing the name and 'mappabile' length of each chromosome.
The mappability tables generated for several species (36 mer and 50 mer) are provided in mptable directory, which are based on the code from [Peakseq](http://info.gersteinlab.org/PeakSeq). See the manual for [DROMPA3](https://github.com/rnakato/DROMPA3) for detail.

# 5. Reference

Nakato R., Shirahige K. Sensitive and robust assessment of ChIP-seq read distribution using a strand-shift profile, Bioinformatics, 2018, https://doi.org/10.1093/bioinformatics/bty137
