#!/usr/bin/env python

"""
Converts .pro format to the same as rnacounter's output.
Doc: http://sammeth.net/confluence/display/SIM/.PRO+Transcriptome+Profile

Usage: parse_pro.py <file.pro> <assembly_name>

"""

from bbcflib.genrep import Assembly
from docopt import docopt

def main(assembly, filename):
    a = Assembly(assembly)
    tmap = a.get_transcript_mapping()

    # Get total num of reads

    f = open(filename)
    g = open('simulation/count_simulation.txt','wb')

    header = ['ID','Count','RPKM','Chrom','Start','End','Strand','GeneName','Length','Type','Sense','Synonyms']
    g.write('\t'.join(header)+'\n')

    for line in f:
        loc, tid, coding, length, \
            expr_fraction, expr_number, lib_fraction, lib_number, seq_fraction, seq_number, \
            cov_fraction, chisq, var_coeff = line.split('\t')
        chrom,coord = loc.split(':')
        start,end = coord[:-1].split('-')
        strand = '1' if coord[-1] == 'W' else '-1'
        nreads = float(seq_number)
        if nreads != 0:
            ntotal = nreads / float(seq_fraction)
            rpkm = 1e9 * nreads / (float(length) * ntotal)
        else:
            rpkm = 0.0
        t = tmap.get(tid)
        if t is not None:
            newline = [tid, seq_number, str(rpkm), chrom,start,end,strand, t.gene_name,length,'transcript','.','.']
            g.write('\t'.join(newline)+'\n')

    f.close()
    g.close()

if __name__ == '__main__':
    args = docopt(__doc__)
    assembly = args['<assembly_name>']
    filename = args['<file.pro>']
    main(assembly, filename)
