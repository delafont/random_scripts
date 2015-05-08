#!/bin/bash

# Count with rnacounter

if [ $# -lt 2 ]
then
  echo "Usage: $0 <bam> <gtf>"
  exit
fi

bam=$1
gtf=$2
superjob='rnac'

#suffix='genes_raw'
#job=${superjob}_${suffix}
#bsub -J $job -e logs/${job}.err -o logs/${job}.out " \
#    rnacounter ${bam} ${gtf} -t genes -m raw > rnacounter/${job}.txt
#"
#
#suffix='genes_indirect-nnls'
#job=${superjob}_${suffix}
#bsub -J $job -e logs/${job}.err -o logs/${job}.out " \
#    rnacounter ${bam} ${gtf} -t genes -m indirect-nnls > rnacounter/${job}.txt
#"

suffix='transcripts_nnls'
job=${superjob}_${suffix}
bsub -J $job -e logs/${job}.err -o logs/${job}.out " \
    rnacounter ${bam} ${gtf} -t transcripts --nh -m nnls > rnacounter/${job}.txt
"
