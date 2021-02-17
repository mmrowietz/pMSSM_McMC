#!/bin/bash

# set up
source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-centos7-gcc8-opt/setup.sh

# copy code
xrdcp root://cmseos.fnal.gov/EOSDIR/pMSSM_McMC.tar.gz .
tar -zxvf pMSSM_McMC.tar.gz
cd pMSSM_McMC

# compile FeynHiggs
cd packages/
rm -rf FeynHiggs-2.16.1 && tar -zxvf FeynHiggs-2.16.1.tar.gz
cd FeynHiggs-2.16.1
./configure && make && make install
cd ../../

# compile SPheno
cd packages/
rm -rf SPheno-4.0.4 && tar -zxvf SPheno-4.0.4.tar.gz
cd SPheno-4.0.4
sed -i "/^F90/c\F90=gfortran" Makefile
make --version && make
cd ../../

# compile superiso
cd packages/
rm -rf superiso_v4.0 && tar -zxvf superiso_v4.0.tgz
cd superiso_v4.0
cp ../slha_chi2_reduced.c .
make
make slha
make slha_chi2
make slha_chi2_reduced
cd ../../

# compile higgsbounds
cd packages/
rm -rf higgsbounds && tar -zxvf higgsbounds.tar.gz
cd higgsbounds
mkdir build && cd build
cmake .. -DFeynHiggs_ROOT=FeynHiggs-2.16.1
make
cd ../../../

# compile higgssignals
cd packages/
rm -rf higgssignals && tar -zxvf higgssignals.tar.gz
cd higgssignals
mkdir build && cd build
cmake .. -DFeynHiggs_ROOT=FeynHiggs-2.16.1
make
cd ../../../

# compile GM2Calc
cd packages/
rm -rf GM2Calc-1.7.3 && tar -zxvf v1.7.3.tar.gz
cd GM2Calc-1.7.3
mkdir build && cd build
cmake ..
make
cd ../../../

# compile micromegas
#cd packages/
#rm -rf micromegas_5.2.4 && tar -zxvf micromegas_5.2.4.tgz
#cd micromegas_5.2.4
#gmake
#cd MSSM
#gmake main=main.c
#cd ../../../

# make output directory
mkdir OUTDIR

# run code
python SCRIPTNAME OPTIONS -o OUTDIR

# move output to eos
xrdcp -rf OUTDIR root://cmseos.fnal.gov/EOSDIR/
