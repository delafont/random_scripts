
import os, sys, csv

def readcel(filename, skip=0):
    f = open(filename, "rb")
    extension = filename.split('.')[-1]
    basename = filename[:-len(extension)]
    g = open(basename+"bed", "wb")
    r = csv.reader(f, delimiter="\t")
    w = csv.writer(g, delimiter="\t")
    i = 0
    g.write("track type=bedGraph name='"+basename+"' description='"+basename+"' \n")
    for line in r:
        if len(line)==2:
            try:
                l = ["chr2", line[0], str(int(line[0])+25), "probe_"+str(i), line[1]]
                w.writerow(l)
                i+=1
            except: pass
    f.close()
    g.close()
    return basename+"bed"
