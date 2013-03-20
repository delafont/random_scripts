# Sort and index all bam files in the current folder
for file in *.bam; do `samtools sort $file $file.sorted`; done  # sort all bam files
ls -I *sorted* | grep bam | xargs rm -r                       # remove all bam files that dont contain 'sorted'
rename '.sorted.bam' '' *                                       # remove all occurences of 'sorted.bam' in all file names
for file in *; do `samtools index $file`; done                  # index all bam files
