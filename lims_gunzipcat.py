#!/usr/bin/env python

''' Small script to concat LIMS data using cat. Note the script for now only works
on single end data.
By Julien Duc, Evarist Planet'''

from subprocess import call
from glob import glob
import sys
import re
import os

if len(sys.argv) > 1 and sys.argv[1] in ['-h','--help']:
    print '''
    Usage: lims_gunzipcat.py
        Will unzip and cat all files ending in .fastq.fz in current directory
        to files with sampleName.fq.
        Using bsub and low memory - This is the best <3

    Args:
        opt:
            --paired: In case of having paired end data.
'''
    sys.exit()

if not os.path.exists('../../logs'):
    os.makedirs('../../logs')

with open("../../logs/lims_gunzipcat.log", "w") as log:
    if '--paired' in sys.argv:
        files = glob('*_R1_*.fastq.gz')
        sample = set([f[0:-31] for f in files])
        for s in sample:
            snew = s[:-1]
            m = re.findall("({snew})_(\w{{12}})_(.*R1.*)\.(fastq.gz)".format(snew=snew), "\n".join(files))
            sinfo = [i[2] for i in m]
            idx = [tmp[0] for tmp in sorted(enumerate(sinfo),key=lambda i:i[1])]
            mo = [m[n][0] + "_" + m[n][1] + "_" + m[n][2] + "." + m[n][3] for n in idx]
            cmd = "zcat " + " ".join(mo) + " > {snew}_R1.fq"
            cmd = cmd.format(snew=snew)
            # print cmd
            log.write(cmd + "\n\n")
            call('bsub "' + cmd + '"', shell=True)
        files = glob('*_R2_*.fastq.gz')
        sample = set([f[0:-31] for f in files])
        for s in sample:
            snew = s[:-1]
            m = re.findall("({snew})_(\w{{12}})_(.*R2.*)\.(fastq.gz)".format(snew=snew), "\n".join(files))
            sinfo = [i[2] for i in m]
            idx = [tmp[0] for tmp in sorted(enumerate(sinfo),key=lambda i:i[1])]
            mo = [m[n][0] + "_" + m[n][1] + "_" + m[n][2] + "." + m[n][3] for n in idx]
            cmd = "zcat " + " ".join(mo) + " > {snew}_R2.fq"
            cmd = cmd.format(snew=snew)
            # print cmd
            log.write(cmd + "\n\n")
            call('bsub "' + cmd + '"', shell=True)
    else:
        files = glob('*.fastq.gz')
        sample = set([f[0:-31] for f in files])
        for s in sample:
            snew = s[:-1]
            m = re.findall("({snew})_(\w{{12}})_(.*R1.*)\.(fastq.gz)".format(snew=snew), "\n".join(files))
            sinfo = [i[2] for i in m]
            idx = [tmp[0] for tmp in sorted(enumerate(sinfo),key=lambda i:i[1])]
            mo = [m[n][0] + "_" + m[n][1] + "_" + m[n][2] + "." + m[n][3] for n in idx]
            cmd = "zcat " + " ".join(mo) + " > {snew}_R1.fq"
            cmd = cmd.format(snew=snew)
            # print cmd
            log.write(cmd + "\n\n")
            call('bsub "' + cmd + '"', shell=True)

