#!/bin/bash

# Count with sailfish

if [ $# -lt 2 ]
then
  echo "Usage: $0 <fastq> <fasta_trx>"
  exit
fi

fastq=$1
fasta_trx=$2

# Load Sailfish in path
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/bin/Sailfish-0.6.3-Linux_x86-64/lib
PATH=~/bin/Sailfish-0.6.3-Linux_x86-64/bin:${PATH}

job='sail'
bsub -q long -J $job -e logs/${job}.err -o logs/${job}.out \
-M 15000000 -R "rusage[mem=15000] span[ptile=16]" "
    if [ ! -d sailfish/index ]; then
        sailfish index --force --threads 16 --transcripts $fasta_trx --kmerSize 20 --out sailfish/index
    fi
    sailfish quant -i sailfish/index/ -p 16 -l 'T=SE:S=U' -r $fastq -o sailfish/sailfish_result
"

