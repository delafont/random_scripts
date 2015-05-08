#!/bin/bash

# Count with cuffquant

if [ $# -lt 2 ]
then
  echo "Usage: $0 <bam> <gtf>"
  exit
fi

bam=$1
gtf=$2

export PATH=/software/bin:$PATH;
module add UHTS/Assembler/cufflinks/2.2.1;

job='cuff'
bsub -J $job -e logs/${job}.err -o logs/${job}.out -M 8000000 -R "rusage[mem=8000] span[ptile=16]" " \
    cuffquant $gtf $bam -o cuffquant/
    cuffnorm -o cuffquant/ -p 16 $gtf cuffquant/abundances.cxb cuffquant/abundances.cxb
"

