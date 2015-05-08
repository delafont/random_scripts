
md5=cbcc5aeeb39d29065c6641aafd5ccaa430706008
gtf=/scratch/cluster/monthly/jdelafon/jerome3/tophat/index_hg38/${md5}.gtf
tophat="/software/UHTS/Aligner/tophat/2.0.13/bin/tophat"

for fastq in fastq/*.fq.gz; do
    group=${fastq%.fq.gz}
    bsub -J $group -e logs/tophat.err -o logs/tophat.out -M 8000000 -R "rusage[mem=8000] span[ptile=6]" " \
        tophat -p 6 -g 1 --no-novel-juncs --no-novel-indels --b2-sensitive -G $gtf  \
        --transcriptome-index tophat/index_hg38/trx  \
        -o tophat/$group \
        tophat/index_hg38/$md5 \
        $fastq
    "
done;

#for group in 1160 1291 68 DelEx GAPO1 GAPO2 Ma San; do
#done
