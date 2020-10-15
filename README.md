# pMSSM_McMC
pMSSM McMC scan code
### description ###
Pure python version of the mcmc code used by me for the CMS pMSSM effort.
An McMC chain can by run by excuting the script mcmc.py with the following options.
The public tools in the package directory need to be compiled.



### non-standard libraries ###
* `argparse`  
* `ROOT`: Database format for point storage  
* `numpy`: to interface python to ROOT  
* `tqdm` (only cosmetic)


### how to run ###
run mcmc.py (python 2.7) with the following options

* `-m, --mode ["new","resume"]`:  
   new: starts a new chain from a random point in the pMSSM, sampling from a flat prior in the 19 pMSSM dimensions.  
    resume: resumes a previously created chain. Specify the root file containing the previous chain. The chain is continued from the point in the input root file that has the highest iteration index. The input file is specified with the `-i, [input]` argument.
* `-i, --input [input]`,If runmode resume, specify input root file (default = None)
* `-o, --output [output]`,Specify an output directory (required)
* `-n, --npoints [npoints]`, How many points to run the MCMC for (default=1000)
* `-c, --chain_index [index]`, chain index for the chain. Points are later identified by the chain_index and the iteration index (default = 1)
* `-mi, --move_interval [interval]` , How many points to generate before starting a new root file and moving the already created points to the output directory (default = 300)



A parallelized batch execution of different chains needs to be implemented and adapted to the local batch system.


### How to add observables ###
1. Add the value for the observable and its uncertainty to the likelihood_contributions dictionary in likelihood.py  
2. Add a way to retrieve the observable to the get_observables function in mcmc.py. This may require you to write an interface to external tools.  
3. If you want to change the way the likelihood is computed for the McMC, check out and change likelihood.py to your needs  


todo:
-mode to start from non-random point
