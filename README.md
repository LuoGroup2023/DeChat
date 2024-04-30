Repeat and haplotype aware error correction in nanopore sequencing reads with DeChat

## Description

Error correction is the canonical first step in long-read sequencing data analysis. Nanopore R10 reads have error rates below 2\%.we introduce DeChat, a novel approach specifically designed for Nanopore R10 reads.DeChat enables repeat- and haplotype-aware error correction, leveraging the strengths of both de Bruijn graphs and variant-aware multiple sequence alignment to create a synergistic approach. This approach avoids read overcorrection, ensuring that variants in repeats and haplotypes are preserved while sequencing errors are accurately corrected.


Dechat is implemented with C++ and Python3.

## Installation and dependencies
DeChat relies on the following dependencies:
- [bcalm]
- [lordec]
- [minimap2]
- [GraphAligner]
- [hifiasm]
- gcc 9.5+ 
- cmake 3.2+
- zlib

#### 1.Install from [Conda]() 
This is easy and recommended:
```
conda create -n dechat
conda activate dechat
conda install -c bioconda dechat
```

#### 2.Install from source code
To install DeChat, firstly, it is recommended to intall the dependencies through [Conda]():
```
conda create -n vechat
conda activate vechat
conda install -c bioconda minimap2 yacrd fpa=0.5
```

## Running and options
If you pulled the source repo; to run dechat 
```bash
git clone https://github.com/LuoGroup2023/DeChat.git
cd dechat
mkdir build &&cd build
cmake ..
make -j 16
```

The input read file is only required and the format should be FASTA/FASTQ (can be compressed with gzip). Other parameters are optional.
Please run `dechat` to get details of optional arguments. 

```
positional arguments:
  sequences             input file in FASTA/FASTQ format (can be compressed
                        with gzip) containing sequences used for correction

optional arguments:
  Input/Output:
    -o STR       prefix of output files [(null)]
    -t INT       number of threads [1]
    -h           show help information
    --version    show version number
  Error correction round 1:
       -r1           number of DBG min [2]
  Error correction round 2:
       -r            second_number_of_round [3]
       -e            max_ov_diff_ec[0]
```

## Examples




