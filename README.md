# SSP (Strand shift profile)

# 1. Overview

# 2. Install
SSP is written in C++ and requires the following programs and libraries:
* [Boost C++ library](http://www.boost.org/)
* [GTK library](http://www.gtk.org/)
* [GNU Scientific Library](http://www.gnu.org/software/gsl/)
* [zlib](http://www.zlib.net/)
* [SAMtools](http://samtools.sourceforge.net/) (for BAM formatted input)

#### 2.1. Install required libraries
for Ubuntu and Debian:

    sudo apt-get install git build-essential libboost-all-dev libgsl-dev libz-dev samtools
 
for CentOS and Red Hat:

    sudo yum -y install git gcc-c++ boost-devel zlib-devel gsl-devel
and install samtools from [the website](http://samtools.sourceforge.net/).

#### 2.3. Install SSP
    git clone https://github.com/rnakato/SSP.git
    cd SSP
    make -j4

#### 2.4. Add the PATH environment variable
For example, if you downloaded SSP into the $HOME/my_chipseq_exp directory, type:

    export PATH = $PATH:$HOME/my_chipseq_exp/SSP/bin

# 3. Usage
The simplest command is:

    ssp -i ChIP.sam -o ChIP --gt genometable.txt
then the output files (prefix: "ChIP") are generated in the directory "sspout (default)".
The format of input file is automatically detected by postfix(.sam/.bam/.bowtie/.tagalign(.gz)). If the detection does not work well, supply -f option (e.g., "-f BAM").

If the mappable genome table is supplied

     ssp -i ChIP.bam -o ChIP --gt genometable.txt --mptable mptable.txt -p 4
"-p 4" specifies the number of CPUs used. Mappable genome table files for several species and mapping parameters are included in data directory.


Multiple input files are allowed (separated by ",")

     ssp -i ChIP1.bam,ChIP2.bam,ChIP3.bam -o ChIP --gt genometable.txt 


Note that the chromosome length should be enough longer than the background length specified. For small genomes (e.g., yeast), the background length should be shorten:

     ssp -i ChIP1.bam -o ChIP --gt genometable.txt --ng_from 10000 --ng_to 50000 --ng_step 500
     
In this parameter, the background region for FCS is the average ranging from 10k to 50k at steps of 500 bp.

In default, FCS is calcutated for 10M reads. If the number of nonredundant reads are smaller than 10M, specify smaller number for fair comparison among samples as follows:

     ssp -i ChIP1.bam -o ChIP --gt genometable.txt --num4ssp 5000000

When specifying smaller read number for --num4ssp, FCS score becomes smaller, but the magnitude relation among samples is consistent.

# 4. Annotation data

# 5. Reference

Under preparation.
