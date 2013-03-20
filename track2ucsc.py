from bbcflib import genrep
import re

def translate_chr_names(file, assembly_id, output=None):
    """
    Bed or bam files as known by Ensembl for instance use chromosome identifiers that
    are not recognized by the UCSC visualizer. This function uses GenRep to translate
    IDs of the form 'NC_312432' to 'chr8'.
    * *file*: input file name.
    * *assembly_id*: GenRep identifier - e.g. 'hg19' or 76 for human.
    * *output*: output file name.
    * *header*: If a header is present, *True* will tell the script to skip the first
    *header_length* lines.
    """
    if not output: output = file+'.chr_renamed'
    o = open(output,'w')
    f = open(file,'r')
    
    g = genrep.GenRep(url='http://bbcftools.vital-it.ch/genrep/',root='/db/genrep')
    assembly = g.assembly(assembly_id)
    chromosomes = assembly.chromosomes

    for line in f:
        found = re.search('[0-9]*_NC_[0-9]*\.[0-9]*',line)
        if found:
            whole = found.group(0)
            chr_id = (int(whole.split('_')[0]),
                      unicode(whole.split('_')[1]+'_'+whole.split('_')[2].split('.')[0]),
                      int(whole.split('.')[1]))
            chr_name = chromosomes[chr_id]['name']
            line = line.replace(whole,chr_name)
        o.write(line)
    f.close()
    o.close()


def add_names(file, output=None):
    if not output: output = file+'.with_names'
    o = open(output,'w')
    f = open(file,'r')

    k=0
    for line in f:
        k+=1
        line = line.strip('\n')+'\t'+'JUNC_'+str(k)+'\n'
        o.write(line)
