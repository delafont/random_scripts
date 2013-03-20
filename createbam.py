
#with open("/archive/epfl/bbcf/jdelafon/mappings/bb27b89826b88823423282438077cdb836e1e6e5.pickle",'rb') as f:
#    import cPickle
#    mappings = cPickle.load(f)
#    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = mappings
#
#gene_name = "ubl4b"
#gene_id = "ENSG00000186150"
#chr = 1
#start = 110655062
#end = 110656569
#strand = 1
#trans = trans_in_gene[gene_id][0] # ENST00000334179
#exon = exons_in_trans[trans][0] # ENSE00001331233
#
#
#gene_name = "cux2"
#gene_id = "ENSG00000111249"
#chr = 12
#start = 111471828
#end = 111788358
#strand = 1
#trans = trans_in_gene[gene_id]# ['ENST00000261726', 'ENST00000397643', 'ENST00000551604', 'ENST00000552889']
#
#gene_name = "1700090g07rik" # not found in ensembl
#
#gene_name = "batf2"
#gene_id = "ENSG00000168062"
#chr = 11
#start = 64755415
#end = 64764517
#strand = -1
#trans = trans_in_gene[gene_id] # ['ENST00000301887','ENST00000435842','ENST00000527454','ENST00000527716','ENST00000534177']
#trans = trans_in_gene[gene_id][0] # ENST00000301887
#exon = exons_in_trans[trans] # ['ENSE00001408542', 'ENSE00001414327', 'ENSE00002151448']
#exon = exons_in_trans[trans][0] # ENSE00001408542
#
#import pysam
#sam = pysam.Samfile("/scratch/cluster/monthly/jdelafon-el/rnaseq_MEF.files/IRYC01kOdpJFIwE0qhJN","rb")
#
#class Counter:
#    cCounts = 0
#    def __call__(self, alignment):
#        self.counts += 1
#
#c = Counter()
#sam.fetch(exon, start,end,callback=c)
#
#
## control
#gene_name = "Gapdh"
#gene_id = "ENSG00000111640"
#chr = 12
#start = 6643093
#end = 6647537
#strand = 1
#trans = trans_in_gene[gene_id]
#"""['ENST00000229239',
#    'ENST00000496049',
#    'ENST00000396856',
#    'ENST00000492719',
#    'ENST00000396861',
#    'ENST00000474249',
#    'ENST00000466588',
#    'ENST00000396859',
#    'ENST00000466525',
#    'ENST00000460556',
#    'ENST00000396858',
#    'ENST00000450282']"""
#exons = []
#for t in trans:
#    exons.extend(exons_in_trans[t])
#
#import pysam
#sam = pysam.Samfile("/scratch/cluster/monthly/jdelafon-el/rnaseq_MEF.files/IRYC01kOdpJFIwE0qhJN","rb")
#
#class Counter:
#    cCounts = 0
#    def __call__(self, alignment):
#        self.counts += 1
#
#counts = {}
#totalcount = 0
#for e in exons[1:]:
#    c = Counter()
#    sam.fetch(e, start,end,callback=c)
#    totalcount += c.counts
#    counts[e] = c.counts
#
#
##--------------
#
#import pysam
#sam = pysam.Samfile("/archive/epfl/bbcf/jdelafon/lims/RNAseq_full.files/ePF0QuyPZ7qrrudG7dqX","rb")
#
## control
#gene_name = "Gapdh"
#gene_id = "ENSG00000111640"
#start = 6643093
#end = 6647537
#exons = []
#for t in trans_in_gene[gene_id]:
#    exons.extend(exons_in_trans[t])
#
#labels = [t['SN'] for t in sam.header['SQ']]
#labels = [l.split("|")[0] for l in labels]
#common = [e for e in exons if e in labels]
#
#def c(alignement):
#    print alignement
#
#for e in exons:
#    sam.fetch(e,start,end,callback=c)


#----------------------

fakebam_content = [
"@HD\tVN:1.0\tSO:unsorted\n",
"@SQ\tSN:ENSE00002188685|ENSG00000111640|6645660|6645759|1\tLN:100\n",
"@SQ\tSN:ENSE00001902446|ENSG00000111640|6647267|6647537|1\tLN:271\n",
"C3PO_0038:5:79:17771:14659#0/1\t16\tENSE00001902446|ENSG00000111640|6647267|6647537|1\t167\t255\t75M\t*\t0\t0\tAATCTCCCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGT\t`]Z`V``bb\edagdedbdbbedbacdb]abbb_bQ]`ZZVcfcc]ffffggcdgggggggggeggggggggggg\tXA:i:0\tMD:Z:75\tNM:i:0\tNH:i:5\n",
"C3PO_0038:5:93:2476:20366#0/1\t16\tENSE00001902446|ENSG00000111640|6647267|6647537|1\t167\t255\t75M\t*\t0\t0\tAATCTCCCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGT\tBBBBBBBBBB__\_V_^^_a\\\XR_aaaYabaaa_aa]_Xchhbaeffdgdfghhhhhghhhhhhhhhhhhhhh\tXA:i:0\tMD:Z:75\tNM:i:0\tNH:i:5\n",
"C3PO_0038:5:72:16306:11928#0/1\t16\tENSE00001902446|ENSG00000111640|6647267|6647537|1\t167\t255\t75M\t*\t0\t0\tAATCTCCCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGT	\[PZP\\`^\^^R^Uaacccccaac__b`Qbbb^b``b``]cffc]fffagddaggggggggggggggggggggg\tXA:i:0\tMD:Z:75\tNM:i:0\tNH:i:5\n",
"C3PO_0038:5:94:4758:4564#0/1\t0\tENSE00002188685|ENSG00000111640|6645660|6645759|1\t6\t255\t75M\t*\t0\t0\tGTCGTATTGGGCGCCTGGTCACCAGGGCTGCTTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCT\thhhghhhhhhhhhhhhhhhhhhhghhhghhghghhgfhhhhhhfhhhhchhdegeggdgghhdghehffcdhffb\tXA:i:0\tMD:Z:75\tNM:i:0\tNH:i:3\n",
"C3PO_0038:5:94:14981:4668#0/1\t0\tENSE00002188685|ENSG00000111640|6645660|6645759|1\t6\t255\t75M\t*\t0\t0\tGTCGTATTGGGCGCCTGGTCACCAGGGCTGCTTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCT\tghhhhhhhfgdhhhhgghfehhhcffhahfhchhhhfhhchchdhhhfahfdRffhabegacWaacc`cadcggd\tXA:i:0\tMD:Z:75\tNM:i:0\tNH:i:3\n" ]

samfile = open("fakebam","w")
samfile.writelines(fakebam_content)
samfile.close()

import pysam
sam = pysam.Samfile("fakebam","r")

class Counter(object):
    def __init__(self):
        self.n = 0
    def __call__(self, alignment):
        self.n += 1

c = Counter()
sam.fetch("ENSE00002188685|ENSG00000111640|6645660|6645759|1",0,100,callback=c)
count = c.n

sam.close()
