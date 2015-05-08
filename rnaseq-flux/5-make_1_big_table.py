#!/usr/bin/env python

## Add annotation
#from bbcflib.genrep import Assembly()
#a = Assembly('hg19')
#tmap = a.get_transcript_mapping()

# Get counts from FluxSimulator
sim_rpkm = {}
sim_count = {}
with open("simulation/count_simulation.txt") as f:
    f.readline()
    for line in f:
        L = line.strip().split('\t')
        enst,count,rpkm,chrom,start,end,strand,gene_name,length,type,sense,synonyms = L
        sim_count[enst] = count
        sim_rpkm[enst] = rpkm

# Get counts from rnacounter and record annotation
rna_rpkm = {}
rna_count = {}
annot = {}
with open("rnacounter/rnac_transcripts_nnls.txt") as f:
    f.readline()
    for line in f:
        L = line.strip().split('\t')
        enst,count,rpkm,chrom,start,end,strand,gene_name,length,type,sense,synonyms = L
        rna_count[enst] = count
        rna_rpkm[enst] = rpkm
        annot[enst] = (chrom,start,end,strand,length,gene_name)
        if synonyms != '.':
            for syn in synonyms.split(','):
                rna_count[syn] = count  # duplicate!
                rna_rpkm[syn] = rpkm    # duplicate!
                annot[syn] = (chrom,start,end,strand,length,gene_name)

# Get counts from sailfish
sail_rpkm = {}
sail_count = {}
if 1:
    with open("sailfish/sailfish_result/quant_bias_corrected.sf") as f:
        [f.readline() for _ in range(5)]
        for line in f:
            enst,length,tpm,rpkm,kpkm,kmers,reads = line.strip().split('\t')
            sail_count[enst] = reads
            sail_rpkm[enst] = rpkm

# Get counts from cuffquant
cuff_fpkm = {}
cuff_count = {}
if 1:
    with open("cuffquant/isoforms.count_table") as f:
        f.readline()
        for line in f:
            enst,count = line.strip().split('\t')[:2]
            cuff_count[enst] = str(round(float(count),4))
    with open("cuffquant/isoforms.fpkm_table") as f:
        f.readline()
        for line in f:
            enst,fpkm = line.strip().split('\t')[:2]
            cuff_fpkm[enst] = fpkm

# Get counts from RSEM
rsem_fpkm = {}
rsem_count = {}
if 1:
    with open("rsem/count_rsem.isoforms.results") as f:
        f.readline()
        for line in f:
            enst,gene,length,efflen,expc,tpm,fpkm,isopct = line.strip().split('\t')
            rsem_fpkm[enst] = fpkm
            rsem_count[enst] = expc #str(1e-9 * 183263572 * float(fpkm) * int(length))


#--------------------------#

# Rewrite the table with counts from other methods - transcripts only
table_tcounts = open("all_counts.txt","wb")
table_tcounts.write('\t'.join(['ENST',
                               'SIMcount','NNLScount','CUFFcount','SAILcount','RSEMcount',
                               'SIMrpkm','NNLSrpkm','CUFFfpkm','SAILrpkm','RSEMfpkm',
                               'Chrom','Start','End','Strand','Length','GeneName'])+'\n')

for enst in sim_count.keys():
    simc = sim_count[enst]
    simr = sim_rpkm[enst]
    rnac = rna_count.get(enst,'NA')
    rnar = rna_rpkm.get(enst,'NA')
    cuffc = cuff_count.get(enst,'NA')
    cuffr = cuff_fpkm.get(enst,'NA')
    sailc = sail_count.get(enst,'NA')
    sailr = sail_rpkm.get(enst,'NA')
    rsemc = rsem_count.get(enst,'NA')
    rsemr = rsem_fpkm.get(enst,'NA')
    chrom,start,end,strand,length,gene_name = annot.get(enst,'NA')
    newline = [enst,
               simc,rnac,cuffc,sailc,rsemc,
               simr,rnar,cuffr,sailr,rsemr,
               chrom,start,end,strand,length,gene_name]
    table_tcounts.write('\t'.join(newline)+'\n')
table_tcounts.close()


