from bbcflib.track import track
from bbcflib.gfminer.common import *
from bbcflib.gfminer.stream import concatenate
import itertools
import os,sys

"""junc2bed <filename> <assembly>
Transforms .junc files from SOAPsplice to BED format."""

def main():
    filename = sys.argv[1]
    assembly = sys.argv[2]
    print filename, assembly
    return to_bed(filename,assembly)

def to_bed(filename,assembly):
    t = track(filename,fields=['chr','start','end','strand','score'],chrmeta=assembly,format='txt')
    # Translate chr names
    s = t.read()
    s1 = map_chromosomes(s, t.assembly.chromosomes)
    # Prepare output bed file
    out = track(filename.rstrip('junc')+'bed', fields=['chr','start','end','name','score','strand'])
    out.make_header({'name':filename,'description':filename})
    mode='append'
    # Add junction names
    c = itertools.count()
    s2 = duplicate(s1,'chr','name')
    s3 = apply(s2,'name',lambda x: 'junction'+str(c.next()))
    # Write
    out.write(s3,mode=mode)
    out.close()

def merge_junc_files(trackList,assembly):
    out = track('all.junc',format='txt',fields=['chr','start','end','strand','score'])
    from bbcflib.genrep import Assembly
    a = Assembly(assembly)
    for c in a.chromosomes:
        tl = [track(t,fields=['chr','start','end','strand','score'],format='txt').read(str(c[0])+'_'+c[1]+'.'+str(c[2]))
              for t in trackList]
        #all = concatenate(tl,remove_duplicates=True)
        all = concatenate(tl,group_by=['chr','start','end'],aggregate={'score':lambda x:sum(x)})
        out.write(all,mode='append')

if __name__ == '__main__':
    sys.exit(main())

