# pMSSM_McMC
pMSSM McMC scan code

Pure python version of the mcmc code used by me for the CMS pMSSM effort.
A McMC chain can by run by excuting the script mcmc.py with the following options.
The tools in the package directory need to be compiled.
The following non-standard python libraries are used:


-argparse
-ROOT: Databse format for point storage
-numpy: to interface python to ROOT
-tqdm (only cosmetic)

run mcmc.py (python 2.7) with the following options

A parallelized batch execution of different chains needs to be implemented and adapted to the local batch system.
