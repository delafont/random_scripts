#!/usr/bin/env python

import os,sys
if len(sys.argv) < 3:
    print "Usage: remove_multimapped_reads <bam> <maxhits>"
    sys.exit(1)


from bbcflib.mapseq import remove_duplicate_reads

bam = sys.argv[1]
maxhits = sys.argv[2]
outbam = remove_duplicate_reads(bam, chromosomes=None, maxhits=maxhits, convert=False)
os.rename(outbam, os.path.splitext(bam)[0]+'_nh'+str(maxhits)+'.bam')

