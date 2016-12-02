# SSP (Strand shift profile)

#1. Overview

#2. Install
#### 2.1. Install required libraries
for Ubuntu:

     sudo apt-get install git build-essential libboost-all-dev samtools
 
for CentOS:

     sudo yum -y install zlib-devel gsl-devel gtk2-devel

#### 2.3. Install SSP
    git clone https://github.com/rnakato/SSP.git
    cd SSP
    make
    

#### 2.4. Add the PATH environment variable
For example, if you downloaded SSP into the $HOME/my_chipseq_exp directory, type:

    export PATH = $PATH:$HOME/my_chipseq_exp/SSP/bin

#3. Usage
The simplest command is:

    ssp -i ChIP.sam -o ChIP --gt genometable.txt
then the output files (prefix: "ChIP") are generated in the directory "sspout (default)".

If the input format is BAM, and mappable genome table is supplied

     ssp -i ChIP.bam -o ChIP -f BAM --gt genometable.txt --mptable mptable.txt -p 4
"-p 4" specifies the number of CPUs used. 

#4. Reference

Under preparation.
