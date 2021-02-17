#!/usr/bin/python

import argparse
import os

script = 'mcmc.py'

loc_base = os.environ['PWD']
logdir = 'test-condor'
outdir = '/store/user/jennetd/snowmass/'

locdir = logdir
os.system('mkdir -p  %s' %locdir)

nchains = 1

for ch in range(1,nchains+1):
 
    prefix = 'output'
    print('Submitting '+prefix)

    condor_templ_file = open(loc_base+"/submit.templ.condor")
    sh_templ_file    = open(loc_base+"/submit.templ.sh")
    
    localcondor = locdir+'/'+prefix+".condor"
    condor_file = open(localcondor,"w")
    for line in condor_templ_file:
        line=line.replace('DIRECTORY',locdir)
        line=line.replace('PREFIX',prefix)
        condor_file.write(line)
    condor_file.close()
    
    localsh=locdir+'/'+prefix+".sh"
    sh_file = open(localsh,"w")
    for line in sh_templ_file:
        line=line.replace('SCRIPTNAME',script)
        line=line.replace('OPTIONS','-m new')
        line=line.replace('OUTDIR',prefix)
        line=line.replace('EOSDIR',outdir)
        sh_file.write(line)
    sh_file.close()

    os.system('chmod u+x '+locdir+'/'+prefix+'.sh')
    
    if (os.path.exists('%s.log'  % localcondor)):
        os.system('rm %s.log' % localcondor)
    os.system('condor_submit %s' % localcondor)

    condor_templ_file.close()
    sh_templ_file.close()

print(nchains,"jobs submitted.")

