#!/bin/bash

# By Julien Duc, Evarist Planet

#BSUB -L /bin/bash
#BSUB -J tophat[1-24]
#BSUB -o logs/tophat_%I.out
#BSUB -e logs/tophat_%I.err
#BSUB -u julien.delafontaine@epfl.ch
#BSUB -N
#BSUB -M 8000000
#BSUB -R rusage[mem=8000]
#BSUB -n 4
#BSUB -R "span[hosts=1]"

module add UHTS/Analysis/samtools/0.1.19;
#module add UHTS/Aligner/bowtie/0.12.9
module add UHTS/Aligner/bowtie2/2.2.1;
module add UHTS/Aligner/tophat/2.0.13;

set -e
## for who is the analysis
who="jerome3"
## Folder name for experiment
exp="tophat"
## Monthly or weekly
when="monthly"
## Where are the data
inputdir="."
[ ! -d $inputdir ] && echo "You are not in the right directory you dumbass" && exit 1
tmp=("dummy")  # because LSF arrays are 0-based or shit like this
fastq=($(ls -1 ${inputdir}/fastq/*.fq.gz))
data=("${tmp[@]}" "${fastq[@]}")
    #if you wanna compute only one
    #data=("empty" "theguy")
sample=`basename ${data[$LSB_JOBINDEX]}`
sample=${sample%.fq}
outputdir="/scratch/cluster/$when/$USER/$who/$exp/${sample}"
outdata=${data[$LSB_JOBINDEX]%.fq}

#clean previous logs
mkdir -p logs
rm -f logs/tophat_${LSB_JOBINDEX}.*

rm -rf $outputdir
mkdir -p $outputdir

org="hg38"
md5="cbcc5aeeb39d29065c6641aafd5ccaa430706008"
reads="/scratch/cluster/$when/$USER/$who/fastq/${data[$LSB_JOBINDEX]}"
index="/scratch/cluster/$when/$USER/$who/$exp/index_$org/$md5"
gtf="/scratch/cluster/$when/$USER/$who/$exp/index_$org/${md5}_ENSEMBL.gtf"

## create an alias for /scratch/local directory
localdir="/scratch/local/daily/${USER}/${LSB_JOBID}_${LSB_JOBINDEX}"
mkdir -p $localdir
cd $localdir

## Actual mapping
echo " >>> Working with ${data[$LSB_JOBINDEX]} <<< "
cmd="tophat -p 4 -g 1 --no-novel-juncs --no-novel-indels -G $gtf --transcriptome-index $exp/index_$org/trx --b2-sensitive -o $localdir $index $reads"
eval $cmd

## Move and clean localdir
cp -rv ./* $outputdir
rm -rfv $localdir



