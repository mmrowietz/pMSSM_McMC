#!/bin/bash
npoints=$1

echo "Sourcing setup scripts"
source /root-6.12.06-build/bin/thisroot.sh

ln -s /afs/cern.ch/work/j/jdickins/snowmass/pMSSM_McMC/packages .
ln -s /afs/cern.ch/work/j/jdickins/snowmass/pMSSM_McMC/*.py .

outdir=output_${npoints}
mkdir ${outdir}
python mcmc.py -m new -o ${outdir} -n ${npoints}
