#!/bin/bash

# Replaces identifiers inside the fasta files to match the 'chr1' format.

if [ $# -lt 1 ]
then
  echo "Usage: $0 <assembly_name>"
  exit
fi

assembly=$1
for file in ./*.fa; do
    bsub -e err.txt -o out.txt "rename_fasta_chr $file $assembly"
done;

