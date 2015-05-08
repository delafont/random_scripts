#!/usr/bin/env python

import sys
if len(sys.argv) < 2:
    print "Usage: header_translation <assembly_name>"
    sys.exit(1)

from bbcflib.genrep import Assembly
assembly = sys.argv[1]
a = Assembly(assembly)

ac2name = {}
for k,v in a.chrmeta.items():
    ac2name[v['ac']] = k

f = open("header.sam")
#g = open("reheader.txt", "wb")
h = open("reheader.sam", "wb")

for line in f:
    L = line.split('\t')
    chrom = L[1].split(':')[1]
    length = L[2].split(':')[1]
    newchrom = ac2name[chrom]
    #g.write('%s\t%s' % (newchrom,length))
    h.write(line.replace(chrom,newchrom))

f.close()
g.close()


