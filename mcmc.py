#this file should contain the mcmc decision function and the point generator
import random
from math import sqrt
import argparse
import os
from ROOT import *
import numpy as np
from tqdm import tqdm
import utils
import likelihood
#set up the parameter ranges
parameter_ranges ={}
for parameter in ["mu","M1","M2"]:
    parameter_ranges[parameter] = (-4000,4000)
for parameter in ["M3","Mq1","Mq3","Mu1","Mu3","Md1","Md3"]:
    parameter_ranges[parameter] = (0,10000)
for parameter in ["Ml1","Mr1","Ml3","Mr3","Mh3",]:
    parameter_ranges[parameter] = (0,4000)
for parameter in ["At","Ab","Al"]:
    parameter_ranges[parameter] = (-7000,7000)
parameter_ranges["tb"] = (2,60)

width_coefficient = 0.1 #the width coefficient of the gaussian for the mcmc step. This coefficient is multiplied by the parameter range to give the width of the gaussian.

homedir = "/nfs/dust/cms/user/mrowietm/python_scan/pMSSM_McMC/"
packagedir = homedir+"packages/"
spnexe = packagedir+"SPheno-4.0.4/bin/SPheno"
fhexe = packagedir+"FeynHiggs-2.16.1/bin/FeynHiggs"
devnull = '>& /dev/null'
#devnull = ""

#containers for the tree branches. Numpy arrays are used as an interface between python types and the root branches
tree_branches = {}
tree_branches["likelihood"]= np.zeros(1,dtype = float)
tree_branches["iteration_index"] = np.zeros(1,dtype = int)
tree_branches["accepted_index"] = np.zeros(1,dtype = int)
tree_branches["chain_index"] = np.zeros(1,dtype = int)
for obs in likelihood.likelihood_contributions:
    tree_branches[obs] = np.zeros(1,dtype=float)
for param in parameter_ranges.keys():
    tree_branches[param] = np.zeros(1,dtype=float)

#sign permutations
def get_sign(signchoice):
    sign_permutation = (0,0,0)#sign_permutation = (sign(mu),sign(M1),sign(M2))
    if signchoice % 8 ==0:
        return {"mu":-1,"M1":-1,"M2":-1}
    if signchoice % 8 ==1:
        return {"mu":1,"M1":1,"M2":1}
    if signchoice % 8 ==2:
        return {"mu":1,"M1":1,"M2":-1}
    if signchoice % 8 ==3:
        return {"mu":1,"M1":-1,"M2":1}
    if signchoice % 8 ==4:
        return {"mu":1,"M1":-1,"M2":-1}
    if signchoice % 8 ==5:
        return {"mu":-1,"M1":1,"M2":1}
    if signchoice % 8 ==6:
        return {"mu":-1,"M1":1,"M2":-1}
    if signchoice % 8 ==7:
        return {"mu":-1,"M1":-1,"M2":1}
    
def generate_point(input_point = {},signchoice = 0):
    # mode 1: generate point from flat prior
    from_scratch = len(input_point) == 0
    signs = get_sign(signchoice)
    output_point = {}
    if from_scratch:
        for parameter,parameterrange in parameter_ranges.items():
            if parameter in signs.keys():
                parametervalue = signs[parameter]*abs(random.uniform(parameterrange[0],parameterrange[1]))
            else:
                parametervalue = random.uniform(parameterrange[0],parameterrange[1])                    
            output_point[parameter] = parametervalue
        output_point["mtop"] = 173.1
        output_point["mbottom"] = 4.18
        output_point["alpha_s"] = 0.1181
    else:
        for parameter,parameterrange in parameter_ranges.items():
            in_range = False
            if parameterrange[0]<0:
                width = 0.5*width_coefficient*(parameterrange[1]-parameterrange[0])
            else:
                width = width_coefficient*(parameterrange[1]-parameterrange[0])
            while not in_range:
                parametervalue = random.gauss(input_point[parameter],width)
                in_range = parametervalue > parameterrange[0] and parametervalue < parameterrange[1]
            output_point[parameter] = parametervalue
        output_point["mtop"]= random.gauss(input_point["mtop"],0.9*width_coefficient)
        output_point["mbottom"]= random.gauss(input_point["mbottom"],((0.18+0.04)/2)*width_coefficient)
        output_point["alpha_s"] = random.gauss(input_point["alpha_s"],0.0011*width_coefficient)
    output_point["scale"] = sqrt(output_point["Mq3"]*output_point["Mu3"])
    return output_point
def run_spheno(inpath):
    cmd = " ".join([spnexe,inpath,devnull])
    os.system(cmd)
    error = open("Messages.out","r").read()
    return len(error) == 0
def run_feynhiggs():
    """
    Replace the SPheno Higgs sector with the FeynHiggs one
    fhin = FeynHiggs output file
    spnin = SPheno output file
    """
    fhin = "SPheno.spc.fh-001"#again, writing a lot of files to disk
    spnin = "SPheno.spc"
    cmd = " ".join([fhexe,spnin,devnull])
    os.system(cmd)# run FeynHiggs
    if not os.path.exists(fhin):#check to see whether FeynHiggs produced an output
        print "could not find Feynhiggs output, skipping point"
        return False
    #fhin part: Get the new Higgs parameters
    masses = {"Mh0":"","MHH":"","MA0":"","MHp":""}#dictionary to collect the FeynHiggs Higgs masses
    dmassblock= "Block"# No replacing needs to be done here, can just append the whole block at the end
    alpha = -1# for the FeynHiggs alpha values
    decaytabs = {"25":"","35":"","36":"","37":""}#dictionary to collect the FeynHiggs decay tables
    with open(fhin,"r") as feynin:
        fhtab = feynin.read()
        blocks = fhtab[:fhtab.find("DECAY")].split("BLOCK")#put all the block in a list
        decays = fhtab[fhtab.find("DECAY"):].split("DECAY")#put the decay tables in a list
        for block in blocks:# loop over the slha blocks
            blocklines = block.split("\n")
            blockname =blocklines[0].strip()
            if blockname == "MASS":#we are only interested in the MASS block
                for bline in blocklines[1:-1]: #iterate over the lines with the particle masses
                    ptc = bline[bline.find("#")+1:].strip()#particle name
                    val = bline.split("   ")[-2].strip()#particle mass
                    if ptc in masses.keys():#only interested in the Higgs bosons
                        masses[ptc] = val.lower()#write value to dictionay
            if blockname == "DMASS":#just append the block
                dmassblock+=block
            if blockname == "ALPHA":#This block only contains the one value
                alpha = blocklines[1]
                alpha = alpha[:alpha.find("#")].strip()
        for decay in decays:#loop over decay blocks for the different particles
            pdgid = decay.split("\n")[0]
            pdgid = pdgid[:pdgid.find(".")-1].strip()#this is the pdgid to which the decay block belongs
            if pdgid in decaytabs:#only interested in Higgs boson decay tables
                decaytabs[pdgid] = "DECAY"+decay#save the respective decay blocks
    #spnin part: find the obsolete Higgs-related parts in the SPheno generated SLHA file
    translate = {"h0":"Mh0","H0":"MHH","A0":"MA0","H)+":"MHp"}# dictionary to translate from SPheno naming convention to FeynHiggs naming convention
    with open(spnin,"r") as tmpin:
        spntab = tmpin.read()#original file content
        newcont = spntab#copy of the original file content. This will be modified with the FeynHiggs masses and decays, then replace the original content
        blocks = spntab[:spntab.find("DECAY")].split("Block")#SPheno SLHA blocks
        decays = spntab[spntab.find("DECAY"):].split("DECAY")#SPheno decay tables
        replaceblock = "\n".join([dmassblock,"# Higgs mixing\n"])#Puts the DMASS block below the MASS block (and above the Higgs mixing part of the SLHA)
        newcont = newcont.replace("# Higgs mixing\n",replaceblock)#insert the DMASS block
        #beware: the order of replacement in this string might matter
        for block in blocks:#loop over SLHA blocks
            blocklines = block.split("\n")
            blockname =blocklines[0][:blocklines[0].find("#")].strip()
            if blockname == "MASS":#only interested in MASS block
                for bline in blocklines[1:-2]:#crop the mass block to relevant part, then loop over it
                    ptc = bline[bline.find("#")+1:].strip()#Particle name of the mass block, beware SPheno naming convention differs from FeynHiggs
                    val = bline.split("  ")[-2].strip()# particle mass
                    if ptc in translate.keys():#Only interested in Higgs boson masses
                        replaceline = bline.replace(val,masses[translate[ptc]])#replace the line containing the mass in SPheno with the FeynHiggs line
                        newcont = newcont.replace(bline,replaceline)#replace the SPheno line
            if blockname == "alpha":#replace the alpha value
                spnalpha = blocklines[1]
                spnalpha = spnalpha[:spnalpha.find("#")].strip()
                replaceline = blocklines[1].replace(spnalpha,alpha)
                newcont = newcont.replace(blocklines[1],replaceline)
        #now replace the decay tables
        for decay in decays:
            pdgid = decay.split("\n")[0]
            pdgid = pdgid[:pdgid.find(".")-1].strip()
            if pdgid in decaytabs:#only interested in Higgs boson tables
                newcont = newcont.replace("DECAY"+decay,decaytabs[pdgid])#replace the decay table for the Higgs bosons
        #write out new file
    with open(spnin,"w") as outfile:
        outfile.write(newcont)#overwrite SPheno file
    os.system("rm "+fhin)# clean up
    return True
def get_observables(slhapath):
    returndict = {}
    with open(slhapath,"r") as slhain:
        slhacont = slhain.read()
        blocks = slhacont.split("Block")[1:]
        #get the masses
        mblock = blocks[14]
        masses = mblock.split("\n")[2:-2]
        returndict["mtop"] = {"value":float(" ".join(masses[0].split()).split()[1])}
        returndict["mhiggs"] = {"value":float(" ".join(masses[4].split()).split()[1])}
        #get the higgs mass uncertainty
        dmblock = blocks[15]
        dmasses = dmblock.split("\n")[1:-1]
        try:
            returndict["mhiggs"]["uncertainty"]=float(" ".join(dmasses[0].split()).split()[1])
        except:
            print " ".join(dmasses[0].split()).split()
            exit()
        #get alpha_s
        smblock = blocks[5]
        sminputs = smblock.split("\n")[1:-1]
        returndict["alpha_s"] = {"value":float(" ".join(sminputs[2].split()).split()[1])}
        returndict["mbottom"] = {"value":float(" ".join(sminputs[4].split()).split()[1])}
    return returndict
def setup_tree(outtree):
    outtree.Branch("slha_file",TString())
    outtree.Branch("iteration_index",tree_branches["iteration_index"],"iteration_index/I")
    outtree.Branch("accepted_index",tree_branches["accepted_index"],"accepted_index/I")
    outtree.Branch("chain_index",tree_branches["chain_index"],"chain_index/I")
    outtree.Branch("likelihood",tree_branches["likelihood"],"likelihood/D")
    for obs in likelihood.likelihood_contributions.keys():
        outtree.Branch(obs,tree_branches[obs],obs+"/D")
    for param in parameter_ranges.keys():
        outtree.Branch(param,tree_branches[param],param+"/D")
        

def prepare_fill(point,outtree):
    point_info = {}
    for key,val in point.items():
        point_info[key] = val


    with open("SPheno.spc","r") as spnin:
        slhacont = spnin.read()
        slhafile = slhacont
        blocks = slhacont.split("Block")[1:]
        params = blocks[4].split("\n")[2:-1]
        point_info["M1"] = float(" ".join(params[0].split()).split()[1])
        point_info["M2"] = float(" ".join(params[1].split()).split()[1])
        point_info["M3"] = float(" ".join(params[2].split()).split()[1])
        point_info["At"] = float(" ".join(params[3].split()).split()[1])
        point_info["Ab"] = float(" ".join(params[4].split()).split()[1])
        point_info["Al"] = float(" ".join(params[5].split()).split()[1])
        point_info["mu"] = float(" ".join(params[6].split()).split()[1])
        point_info["tb"] = float(" ".join(params[7].split()).split()[1])
        point_info["Mh3"] = float(" ".join(params[8].split()).split()[1])
        point_info["Ml1"] = float(" ".join(params[9].split()).split()[1])
        point_info["Ml2"] = float(" ".join(params[10].split()).split()[1])
        point_info["Ml3"] = float(" ".join(params[11].split()).split()[1])
        point_info["Mr1"] = float(" ".join(params[12].split()).split()[1])
        point_info["Mr2"] = float(" ".join(params[13].split()).split()[1])
        point_info["Mr3"] = float(" ".join(params[14].split()).split()[1])
        point_info["Mq1"] = float(" ".join(params[15].split()).split()[1])
        point_info["Mq2"] = float(" ".join(params[16].split()).split()[1])
        point_info["Mq3"] = float(" ".join(params[17].split()).split()[1])
        point_info["Mu1"] = float(" ".join(params[18].split()).split()[1])
        point_info["Mu2"] = float(" ".join(params[19].split()).split()[1])
        point_info["Mu3"] = float(" ".join(params[20].split()).split()[1])
        point_info["Md1"] = float(" ".join(params[21].split()).split()[1])
        point_info["Md2"] = float(" ".join(params[22].split()).split()[1])
        point_info["Md3"] = float(" ".join(params[23].split()).split()[1])
        sms =  blocks[5].split("\n")[1:-1]
        point_info["alpha_s"] = float(" ".join(sms[2].split()).split()[1])
        point_info["mbottom"] = float(" ".join(sms[4].split()).split()[1])
        point_info["mtop"] = float(" ".join(sms[5].split()).split()[1])
        masses = blocks[14].split("\n")[2:-2]
        point_info["mhiggs"] = float(" ".join(masses[4].split()).split()[1])
        slha_file = TString(slhafile)
        outtree.SetBranchAddress("slha_file",slha_file)
    for key in tree_branches.keys():
        tree_branches[key][0]=point_info[key]
    return point_info
def run(arguments):
    mode = arguments.mode
    chainix = arguments.chain_index
    tend = arguments.npoints
    start = 0
    move_every = arguments.move_interval
    outdir = arguments.output
    inpath = arguments.input
    lastaccepted = {}
    #select mode, use argparse for this
    move = 1
    #run the thing
    #mode 1: start from new point
    if mode == "new":
        outname = "pMSSM_MCMC_"+str(chainix)+"_"+str(start)+"to"+str(min(start+move_every,start+tend))+".root"
        outroot = TFile(outname,"recreate")
        outtree = TTree("mcmc","mcmc")
        setup_tree(outtree)
        tree_branches["chain_index"][0] = chainix

        finite_lh = False
        signchoice = random.randint(0,7)
        while not finite_lh:
            utils.clean()
            spnerr = False
            while not spnerr:#find a viable point
                utils.clean()
                candidate = generate_point(signchoice = signchoice)#generate a point from flat prior
                spnin = utils.write_spheno_input(candidate)#write the input for spheno
                spnerr = run_spheno(spnin) #run spheno, check if viable point

            if not run_feynhiggs():#run feynhiggs, replace higgs sector
                continue
            observables = get_observables(slhapath = "SPheno.spc") #get observables for the likelihood
            _l = likelihood.get_likelihood(observables)#get likelihood
            finite_lh = _l != 0
        lastaccepted["likelihood"]=_l
        lastaccepted["iteration_index"] = 1
        lastaccepted["accepted_index"] = 1
        lastaccepted["chain_index"] = chainix
        #write point to root, start loop
        for obs in observables.keys():
            lastaccepted[obs] = observables[obs]["value"]
        lastaccepted = prepare_fill(lastaccepted,outtree)#add the rest of the point info, fill the tree branches
        outtree.Fill()
    #mode 2: continue from previous point/root file?
    elif mode == "resume":
        lastaccepted = utils.get_point_from_rootfile(inpath,chainix)
    start = lastaccepted["iteration_index"]+1
    if mode == "resume":
        outname = "pMSSM_MCMC_"+str(chainix)+"_"+str(start)+"to"+str(min(start+move_every,start+tend))+".root"
        print "Creating file "+outname
        outroot = TFile(outname,"recreate")
        outtree = TTree("mcmc","mcmc")
        setup_tree(outtree)
        tree_branches["chain_index"][0] = chainix

    #run
    for iter_ix in tqdm(range(start,start+tend+1)):
        if move == move_every-1 and iter_ix < start+tend-1:
            outtree.BuildIndex("chain_index","iteration_index")
            outtree.Write()
            outroot.Close()
            print "Made "+str(move_every)+" iterations, moving "+outname+" to storage"
            os.system(" ".join(["mv",outname,outdir]))
            outname = "pMSSM_MCMC_"+str(chainix)+"_"+str(iter_ix)+"to"+str(min(iter_ix+move_every,start+tend+1))+".root"
            print "Creating file "+outname
            outroot = TFile(outname,"recreate")
            outtree = TTree("mcmc","mcmc")
            setup_tree(outtree)
            move = -1
        finite_lh = False
        while not finite_lh:
            utils.clean()
            spnerr = False
            while not spnerr:#find a viable point
                utils.clean()
                candidate = generate_point(lastaccepted)#generate a point from the last point
                spnin = utils.write_spheno_input(candidate)#write the input for spheno
                spnerr = run_spheno(spnin) #run spheno, check if viable point
            if not run_feynhiggs():#run feynhiggs, replace higgs sector
                continue
            observables = get_observables(slhapath = "SPheno.spc") #get observables for the likelihood
            _l = likelihood.make_decision(observables,lastaccepted["likelihood"])#get likelihood
            finite_lh = _l != 0
        if _l<0:
            move +=1
            if iter_ix == start+tend:
                print "Made all "+str(tend)+" iterations, moving "+outname+" to storage"
                outtree.BuildIndex("chain_index","iteration_index")
                outtree.Write()
                outroot.Close()
                os.system(" ".join(["mv",outname,outdir]))
            continue #point was not accepted
        lastaccepted["likelihood"]=_l
        lastaccepted["iteration_index"] = iter_ix
        lastaccepted["accepted_index"] =lastaccepted["accepted_index"]+1
        lastaccepted["chain_index"] = chainix
        #write point to root, start loop
        lastaccepted = prepare_fill(lastaccepted,outtree)#add the rest of the point info, fill the tree branches
        outtree.Fill()
        if iter_ix == start+tend:
            print "Made all "+str(tend)+" iterations, moving "+outname+" to storage"
            outtree.BuildIndex("chain_index","iteration_index")
            outtree.Write()
            outroot.Close()
            os.system(" ".join(["mv",outname,outdir]))
        move+=1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--mode",choices=["new","resume"],help="Specify run mode, either new or resume",required=True)
    parser.add_argument("-i","--input",help="If runmode resume, specify input root file",default = None)
    parser.add_argument("-o","--output",help="Specify an output directory",required=True)
    parser.add_argument("-n","--npoints",help = "How many points to run the MCMC for",default=1000,type=int)
    parser.add_argument("-c","--chain_index",help = "chain index for the chain",type = int,default = 1)
    parser.add_argument("-mi","--move_interval",default = 300,help = "How many points to generate before starting a new root file",type=int)
    args=parser.parse_args()
    if args.mode =="resume":
        if args.input == None:
            parser.error("Need to specify an input file in order to resume")
        if args.chain_index!=None:
            print "Ignoring specified chain index, taking from input file"            
    if args.mode == "new":
        if args.input!=None:
            print "Ignoring input file, as we are starting a new chain"
    run(args)

    

