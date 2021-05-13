#!/bin/bash
npoints=$1

echo "Sourcing setup scripts"
source /root-6.12.06-build/bin/thisroot.sh

#git clone https://github.com/jennetd/pMSSM_McMC.git
tar -zxvf pMSSM_McMC-develop.tar.gz
cd pMSSM_McMC

./install-hb.sh 
export HiggsBounds_DIR=${PWD}/packages/higgsbounds/build/

./install-hs.sh 

outdir=output_${npoints}
mkdir ${outdir}
python mcmc.py -m new -o ${outdir} -n ${npoints}
