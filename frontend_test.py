
from bbcflib import frontend
import time, datetime


reload(frontend)
# Creation of Job by Frontend
f = frontend.Frontend(url='http://htsstation.vital-it.ch/rnaseq/')
job = f.job('wApOfMlXGSKBuCfHjw73')
job.groups.keys() == [16,17]
job.groups[16].runs.keys() == [29]
job.groups[17].runs.keys() == [30]


reload(frontend)
# Creation of Job by hand
job = frontend.Job({'id':4, 'created_at':time.localtime(), 'key':'k3y', 'assembly_id':'hg19',
             'description':'desc', 'email':'email', 'options':{'discard':'True'}}  )
job = frontend.Job({'id':'atr', 'created_at':time.localtime(), 'key':'k3y', 'assembly_id':'hg19',
             'description':'desc', 'email':'email', 'options':{'discard':'True'}}  )
job.add_group(frontend.Group({'id':33, 'created_at':time.localtime(),
                     'name':'groupA', 'control':True}))
g33 = job.groups[33]
g33.add_run(frontend.Run({'id':44, 'facility':'fac',
                 'facility_location':'facloc','machine':'mach',
                 'machine_id':13,'lane':1,'url':'url','key':'key'}))
r44 = g33.runs[44]
r44.job_id == 4
r44.group_id == 33


reload(frontend)
# Creation of Job by config file
job,gl = frontend.parseConfig('config_files/jobtest.txt')
job.groups.keys() == [12,22]
job.groups[12].runs.keys() == [13]
job.groups[22].runs.keys() == [33]
