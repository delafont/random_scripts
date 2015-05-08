#!/bin/bash

# Count with rsem

#Get rsem
#   wget http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.19.tar.gz
#   tar xzvf rsem-1.2.19.tar.gz
#   rm rsem-1.2.19.tar.gz
#   cd rsem-1.2.19 && make && cd ..

PATH=~/bin/rsem-1.2.19:$PATH

if [ $# -lt 3 ]
then
  echo "Usage: $0 <genome_fasta> <gtf> <fastq>"
  exit
fi

genome_fasta=$1
gtf=$2
fastq=$3

job='rsem'
bsub -q long -J $job -e logs/${job}.err -o logs/${job}.out \
-M 40000000 -R "rusage[mem=40000] span[ptile=16]" "
    if ! ls rsem/ref.* 1> /dev/null 2>&1; then
        rsem-prepare-reference --bowtie --gtf $gtf $genome_fasta rsem/ref
    fi
    rsem-calculate-expression -p 16 $fastq rsem/ref rsem/$job
"
