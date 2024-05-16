Repeat and haplotype aware error correction in nanopore sequencing reads with DeChat

## Description

Error correction is the canonical first step in long-read sequencing data analysis. Nanopore R10 reads have error rates below 2\%. we introduce DeChat, a novel approach specifically designed for Nanopore R10 reads.DeChat enables repeat- and haplotype-aware error correction, leveraging the strengths of both de Bruijn graphs and variant-aware multiple sequence alignment to create a synergistic approach. This approach avoids read overcorrection, ensuring that variants in repeats and haplotypes are preserved while sequencing errors are accurately corrected.

DeChat can use HIFi or NGS to correct ONT now

Dechat is implemented with C++.

## Installation and dependencies
DeChat relies on the following dependencies:
- gcc 9.5+ 
- cmake 3.2+
- zlib
- boost 1.67

#### 1.Install from [Conda]() 
This is easy and recommended:
```
conda create -n dechat
conda activate dechat
conda install -c bioconda dechat
```
#### 2.Install from source code 
```bash
git clone https://github.com/LuoGroup2023/DeChat.git
conda create -n dechat boost=1.67.0
conda activate dechat
```

## Running and options
If you pulled the source repo; to run dechat 
```bash
cd dechat
mkdir build &&cd build
cmake ..
make -j 16
cd ..
./bin/dechat
```

The input read file is only required and the format should be FASTA/FASTQ (can be compressed with gzip). Other parameters are optional.
Please run `dechat` to get details of optional arguments. 

```
Repeat and haplotype aware error correction in nanopore sequencing reads with DeChat
Usage: dechat [options] -o <output> -t <thread>  -i <reads> <...>
Options:
  Input/Output:
       -o STR       prefix of output files [(null)]
                   The output for the stage 1 of correction is "recorrected.fa", 
                    The final corrected file is "file name".ec.fa;
       -t INT       number of threads [1]
       -h           show help information
       -v --version show version number
       -i           input reads file
       -k INT       k-mer length (must be <64) [21]
  Error correction stage 1 (dBG):
       -r1           set the maximal abundance threshold for a k-mer in dBG [2]
       -d           input reads file for building dBG (Default use input ONT reads) 
  Error correction stage 2 (MSA):
       -r            round of correction in alignment [3]
       -e            maximum allowed error rate used for filtering overlaps [0.04]     
```

## Examples

The example folder contains test data, including the 10X depth sim-ont10.4 data of Escherichia coli diploid and its corresponding reference sequence. The running method is as follows:
```
cd example
dechat -i reads.fa.gz -o reads -t 8
```
### Using HIFi or NGS to correct ONT
```
cd example
dechat -i reads.fa.gz -o reads -t 8 -d HiFi-reads.fq.gz/NGS-reads.fq.g
```


