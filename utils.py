from math import exp,pi
from ROOT import *
import os
def gauss(x,mu,sigma):
    return (1./(2*pi*sigma))*exp(-((x-mu)**2)/(2*sigma**2))
def gauss_pm(x,mu,sigma_p,sigma_m):
    if x<mu:
       return (1./(2*pi*sigma_m))*exp(-((x-mu)**2)/(2*sigma_m**2))
    else:
        return (1./(2*pi*sigma_p))*exp(-((x-mu)**2)/(2*sigma_p**2))

def write_spheno_input(candidate_point):
    #write input file for spheno
    block_modsel = ["BLOCK MODSEL #Model selection","1 0 #General MSSM"]
    block_sphenoinput = ["BLOCK SPhenoInput #Spheno specific input"]
    block_sphenoinput.append("1 -1 #error level")
    block_sphenoinput.append("2 0 #SPA convention")
    block_sphenoinput.append("11 1 #calculate branching ratios")
    block_sphenoinput.append("12 1.0E-04 #write only branching ratios larger than this value")
    block_sphenoinput.append("21 0 #calculate cross sections")
    block_sphenoinput.append("38 3 #use three-loop SM-RGEs and two-loop SUSY-RGEs")
    block_sminputs = ["BLOCK SMINPUTS #Standard model inputs"]
    block_sminputs.append("2 1.166379E-5 #G_F")
    block_sminputs.append("3 "+str(candidate_point["alpha_s"])+" #alpha_s(MZ) SM MSbar")
    block_sminputs.append("4 91.1876 #Z pole mass")
    block_sminputs.append("5 "+str(candidate_point["mbottom"])+" #m_b(m_b) SM MSbar")
    block_sminputs.append("6 "+str(candidate_point["mtop"])+" #m_top(pole)")
    block_extpar = ["BLOCK EXTPAR #Input parameters"]
    block_extpar.append("0 "+str(candidate_point["scale"])+" #Input scale")
    block_extpar.append("1 "+str(candidate_point["M1"])+" #M1")
    block_extpar.append("2 "+str(candidate_point["M2"])+" #M2")
    block_extpar.append("3 "+str(candidate_point["M3"])+" #M3")
    block_extpar.append("11 "+str(candidate_point["At"])+" #At")
    block_extpar.append("12 "+str(candidate_point["Ab"])+" #Ab")
    block_extpar.append("13 "+str(candidate_point["Al"])+" #Al")
    block_extpar.append("23 "+str(candidate_point["mu"])+" #mu")
    block_extpar.append("25 "+str(candidate_point["tb"])+" #tan beta")
    block_extpar.append("26 "+str(candidate_point["Mh3"])+" #mA")
    block_extpar.append("31 "+str(candidate_point["Ml1"])+" #Ml1")
    block_extpar.append("32 "+str(candidate_point["Ml1"])+" #Ml2")
    block_extpar.append("33 "+str(candidate_point["Ml3"])+" #Ml3")
    block_extpar.append("34 "+str(candidate_point["Mr1"])+" #Mr1")
    block_extpar.append("35 "+str(candidate_point["Mr1"])+" #Mr2")
    block_extpar.append("36 "+str(candidate_point["Mr3"])+" #Mr3")
    block_extpar.append("41 "+str(candidate_point["Mq1"])+" #Mq1")
    block_extpar.append("42 "+str(candidate_point["Mq1"])+" #Mq2")
    block_extpar.append("43 "+str(candidate_point["Mq3"])+" #Mq3")
    block_extpar.append("44 "+str(candidate_point["Mu1"])+" #Mu1")
    block_extpar.append("45 "+str(candidate_point["Mu1"])+" #Mu2")
    block_extpar.append("46 "+str(candidate_point["Mu3"])+" #Mu3")
    block_extpar.append("47 "+str(candidate_point["Md1"])+" #Md1")
    block_extpar.append("48 "+str(candidate_point["Md1"])+" #Md2")
    block_extpar.append("49 "+str(candidate_point["Md3"])+" #Md3")
    outstring = "\n".join(["\n".join(block_modsel),"\n".join(block_sphenoinput),"\n".join(block_sminputs),"\n".join(block_extpar)])
    outfile = "genpoint.slha"
    with open(outfile,"w") as SpnIn:
        SpnIn.write(outstring)
    return outfile

def clean():
    if os.path.exists("SPheno.spc"):
        os.system("rm SPheno.spc")
    if os.path.exists("SPheno.spc.fh-001"):
        os.system("rm SPheno.spc.fh-001")
    if os.path.exists("Messages.out"):
        os.system("rm Messages.out")


def get_point_from_rootfile(inpath,chainix):
    inroot = TFile(inpath)
    intree = inroot.Get("mcmc")
    maxiter = int(intree.GetMaximum("iteration_index"))
    intree.GetEntryWithIndex(chainix,maxiter)
    returndict = {}
    for branch in intree.GetListOfBranches():
        if branch.GetName()=="slha_file":
            continue
        returndict[branch.GetName()] = getattr(intree,branch.GetName())
    inroot.Close()
    return returndict
