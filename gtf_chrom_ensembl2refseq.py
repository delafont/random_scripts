

from bbcflib.genrep import Assembly
a = Assembly('hg38')
chrmeta = a.chrmeta

md5="cbcc5aeeb39d29065c6641aafd5ccaa430706008"

filename = "%s_ENSEMBL.gtf" % md5
to = "%s_REFSEQ.gtf" % md5
f = open(filename)
g = open(to, "wb")
for line in f:
    L = line.split('\t')
    ensembl = L[0]
    refseq = chrmeta[ensembl]['ac']
    newline = [refseq]+L[1:]
    g.write('\t'.join(newline))
f.close()
g.close()

