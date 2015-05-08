#!/bin/bash

# Convert sam to bam.
# Rename chromosomes in header; remove unmapped; add NH flag; keep only NH=1.

if [ $# -lt 1 ]
then
  echo "Usage: $0 <assembly_name> [name]"
  exit
else
  assembly=$1 ;
    if [ $# -lt 2 ]
    then
        name='realignment'
    else
        name=$2
    fi
fi

archpath='/scratch/cluster/monthly/jdelafon/rnaseq/FluxSimulator'
path1='bowtie'

# Rename chromosomes: prepare sam header
#head -50 ${path1}/${name}.sam | grep 'SN:' > ${path1}/header.sam
#python ${archpath}/header_translation.py ${assembly}

job='sam2bam'
maxhits='5'
bsub -J ${job} -e logs/${job}.err -o logs/${job}.out "
    samtools view -Shb ${path1}/${name}.sam > ${path1}/${name}.bam
    samtools view -hb -F 4 ${path1}/${name}.bam > ${path1}/${name}_mapped.bam
    add_nh_flag ${path1}/${name}_mapped.bam ${path1}/${name}_nh.bam
    python ${archpath}/remove_multimapped_reads.py ${path1}/${name}_nh.bam $maxhits
    samtools sort ${path1}/${name}_nh_nh${maxhits}.bam ${path1}/${name}_final
    samtools index ${path1}/${name}_final.bam
"

   #samtools reheader ${path1}/reheader.sam ${path1}/${name}.bam > ${path1}/${name}_renamed.bam
   #samtools view -hb -F 4 ${path1}/${name}_renamed.bam > ${path1}/${name}_mapped.bam

   #samtools view -hb -q 1 ${path1}/${name}_nh.bam > ${path1}/${name}_unique.bam
    #samtools sort ${path1}/${name}_nh_unique.bam ${path1}/${name}_final
    #samtools index ${path1}/${name}_final.bam


