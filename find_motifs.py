from bbcflib.track import track
from bbcflib.gfminer.common import sorted_stream, select, cobble, apply
import os, shutil


def fimo(motifs,fasta,qval=True):
    # Run Fimo
    if qval:
        options = "--max-stored-scores 1000000 --verbosity 1 --thresh 0.01 --qv-thresh"
    else:
        options = "--max-stored-scores 1000000 --verbosity 1 --thresh 0.000001"
    cmd = "fimo " + options + " %s %s" % (motifs, fasta)
    print "Running >>",cmd
    os.system(cmd)
    os.system("sort -k2,2n -k3,3n -k4,4n fimo_out/fimo.txt > fimo.txt")

    # Bed output
    t = track('fimo.txt', fields=["name","chr","start","end","strand","score","p-value","q-value","sequence"])
    t.fields = ["name","chr","start","end","strand","a","score","q","sequence"]
    s = t.read()
    s = select(s,['chr','start','end','name','score','strand'])
    s = apply(s,'chr',lambda x:x.split('|')[1])
    s = sorted_stream(s)
    s = cobble(s)
    s = apply(s,'name',lambda x:'|'.join(list(set(x.split('|')))))
    outname = 'fimo.bed'
    bed = track(outname,fields=s.fields)
    bed.make_header(name="TSS_motifs", description="Motifs +-XKb around TSS", mode='overwrite')
    bed.write(s)
    if os.path.exists("fimo_out"): shutil.rmtree("fimo_out")

if 0:
    motifs = "pwm/transfac2012.meme"
    fasta = "../Dlgap1_ensembl_longest.fa"
    fimo(motifs,fasta,qval=False)

###############################################################################

def unique_motifs_list(fimo_bed):
    """Print the list of unique motifs found in the file created above"""
    t = track(fimo_bed)
    name_idx = t.fields.index('name')
    motifs = {}
    for x in t.read():
        name = x[name_idx]
        motifs[name] = motifs.get(name,0) + 1
    print motifs.keys()
    for k,v in sorted(motifs.iteritems()):
        print "%s: %s" % (k,v)

if 0:
    unique_motifs_list('fimo.bed')



