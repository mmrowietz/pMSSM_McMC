from ROOT import *
import random
import argparse
import numpy
"""
Sample from the McMC according to the sampling function desired
how it works:
Iterate through every point in the scan. For each point, calculate a probability to pick the point. 
"""
#define the base probability here
p_base = 0.01

#the sampling factor is used to modify the sampling rate
sampling_factor = 1

sampling_contributions = {} #define which regions to sample by which factor
sampling_contributions["example"] = 5
sampling_contributions["example_2"] = 0.1
p_max = p_base*20#if possible, give the maximum probability to pick a point here to speed up loop
#input
infile = TFile(args.infile)
intree = infile.Get("mcmc")
nEntries = intree.GetEntries()
#output
outfile = TFile(args.outfile,"recreate")
outtree = intree.CloneTree(0)#clone the tree structure from the input tree.

#new branches for output tree
tree_branches = {}
tree_branches["p_pick"] = {"container":np.zeros(0,dtype),"dtype":"D"}#track the pick probabilty for each picked points in order to reweight later

for branch in tree_branches.keys():
    if tree_branches[branch]["dtype"] == "TString":#the container in this case is just a placeholder, since a new TString is created for each new tree entry. I am not sure if this can be done differently
        outtree.Branch(branch,tree_branches[branch]["container"])
    else:
        outtree.Branch(branch,tree_branches[branch]["container"],branch+"/"+tree_branches[branch]["dtype"])


    

#aux functions
def apply_positive_modifiers():
    """
    Apply sampling contributions that increase the sampling rate / pick probability
    """
    f_oversampling = 1
    return f_oversampling
def apply_negative_modifiers(threshold,p_pick):
    """
    Apply sampling contributions that decrease the sampling rate / pick probabilty.
    Order contributions in order of increasing computation time, or decreasing impact on sampling rate.
    After each contribution, check whether threshold>p_pick. If not, the point is not picked and the rest of the contribution do not need to be checked

    """
    f_undersampling = 1
    if apply_example_2():
        f_undersampling*= sampling_contributions["example_2"]

    if threshold>p_pick*f_undersampling: #early skip decision possible
        return f_undersampling

    return f_undersampling # return the undersampling factor
    
    
#main loop
for iEntry in range(nEntries):
    intree.Get(iEntry)
    p_pick = p_base
    #sample a random float between 0 and 1. If the pick probability to pick a point exceeds the pick threshold, it is picked, otherwise it is skipped.
    pick_threshold = random.random()
    if pick_threshold > p_max: #before calculating anything, check if the pick threshold exceeds the maximum possible pick probability
        continue
    #first apply all the positive modifiers
    p_pick *= apply_positive_modifiers()
    #p_pick is now as high as it can get for this point. If it does not exceed the threshold, we can skip calculating the rest of the modifiers
    if pick_threshold > p_pick:
        continue
    #now apply the negative modifiers
    p_pick *= apply_negative_modifiers(pick_threshold,p_pick)
    if pick_threshold > p_pick:
        continue
    outtree.Fill()

sampling_fraction = float(outtree.GetEntries())/nEntries
print "My best estimate for the actual sampling fraction is "+str(sampling_fraction)
    

