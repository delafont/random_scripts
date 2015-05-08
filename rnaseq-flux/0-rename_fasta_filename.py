#!/usr/bin/env python

# Renames the fasta files so that each <name>.fa corresponds to the chromosome <name>.

import os,sys
if len(sys.argv) < 2:
    print "Usage: python rename_fasta_filename.py <assembly_name>"
    sys.exit(1)

from bbcflib.genrep import Assembly
assembly = sys.argv[1]
a = Assembly(assembly)

for k,v in a.chrmeta.items():
    accession = v['ac']   # 10620_NC_000067.6
    name = k
    print name, accession
    if os.path.exists(accession+'.fa'):
        os.rename(accession+'.fa', name+'.fa')



