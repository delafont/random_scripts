#!/usr/bin/env python

"""Generate a fake fastq file from genomic sequences in fasta format."""
#bowtie /db/genrep/nr_assemblies/bowtie/bb27b89826b88823423282438077cdb836e1e6e5 testing_files/test.fastq
import os,sys,random

def set_id(type,paired):
    instrument_name = '100R'
    flowcell_lane = '1'
    tile_number = str(random.randint(1,1000))
    x_coord = str(random.randint(1,10000))
    y_coord = str(random.randint(1,10000))
    id = ':'.join(['EAS'+instrument_name,flowcell_lane,tile_number,x_coord,y_coord])
    if type=='illumina':
        id = 'HWUSI-' + id
        if paired: 
            multiplex = '#0'
            pair = '/1'
            id += multiplex + pair
    else:
        if paired:
            fails_filter = 'N'
            control_bits = '18'
            index_seq = 'ATGC'
            pair = '1'
            id += ' ' + ':'.join([pair,fails_filter,control_bits,index_seq])
    return id

def main(type='illumina', paired=False):
    try:
        args = sys.argv
        paired = len(args)==1
        fasta_name = args[1]
        if fasta_name.endswith('.fa'): extension = '.fa'
        elif fasta_name.endswith('.fasta'): extension = '.fasta'
        else: extension = ''
        fastq_name = fasta_name[:-len(extension)]+'.fastq'
        fastq = open(fastq_name,'wb')
        id = ''
        fasta = open(fasta_name,'rb')
        for line in fasta:
            id = set_id(type,paired)
            quality = '~'*(len(line)-1)
            fastq.write('@'+id+'\n')
            fastq.write(line)
            fastq.write('+'+id+'\n')
            fastq.write(quality+'\n')
        fastq.close(); fasta.close()
        print 'Out:', fastq_name
    except Exception, e:
        raise ValueError("\n\n\tUsage: fasta2fastq <sequences.fa>"+'\n\tRaised: '+str(e)+'\n')

if __name__ == '__main__':
    sys.exit(main())

