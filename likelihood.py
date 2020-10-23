import utils
import random
import math


likelihood_contributions ={}
likelihood_contributions["mtop"] = {"value":173.1,"uncertainty":0.9}
likelihood_contributions["mbottom"] = {"value":4.18,"uncertainty":[0.03,0.04]}
likelihood_contributions["alpha_s"] = {"value":0.1181,"uncertainty":0.0011}
likelihood_contributions["mhiggs"] = {"value":125.26}#good practise would be to cite from where the value is taken
#superiso gives only the SUSY contribution to amu, thus it must be compared against delta amu
likelihood_contributions["a_muon"]={"value":26.8E-10,"uncertainty":4.3E-10}#http://pdg.lbl.gov/2018/reviews/rpp2018-rev-g-2-muon-anom-mag-moment.pdf
likelihood_contributions["BR_B_to_tau_nu"]={"value":1.44E-04,"uncertainty":0.31E-04}#http://www.slac.stanford.edu/xorg/hflav/rare/May2018/radll/OUTPUT/TABLES/PDF/radll.pdf
likelihood_contributions["BR_Ds_to_tau_nu"]={"value":5.48E-02,"uncertainty":0.23E-02}#http://pdg.lbl.gov/2018/listings/rpp2018-list-Ds-plus-minus.pdf
likelihood_contributions["BR_Ds_to_mu_nu"]={"value":5.5E-03,"uncertainty":0.23E-03}#http://pdg.lbl.gov/2018/listings/rpp2018-list-Ds-plus-minus.pdf
#take this from SPheno
#likelihood_contributions["drho"]={"value":3.9E-04,"uncertainty":1.9E-04}#http://pdg.lbl.gov/2018/reviews/rpp2018-rev-standard-model.pdf

#add superiso variables here. Maken the keys root compatible, then make an interface for superiso
def get_likelihood(observables):
    """
    get the likelihood for a point.
    @param observables: Dictionary containing the observables for which the likelihood is to be computed. Observable values are keyed by observable name. If a theory uncertainty is given, it is taken to be the only uncertainty pertaining to the observable.
    """
    product_likelihood = 1

    for obs,obsval in observables.items():
        if "special_case" in obsval.keys():#in case you need a non-standard handling of the likelihood, add the flag "special_case" as one of the keys in the likelihood_contributions dictionary
            continue
        if obs not in likelihood_contributions:
            print(str(obs)+" has no experimental result in database")
            continue
        if "uncertainty" in obsval.keys(): #higgs-type case: center gauss on theory value with width= theory uncertainty, evaluate at experimental value
            product_likelihood *= utils.gauss(likelihood_contributions[obs]["value"],obsval["value"],0.5*obsval["uncertainty"])
        else:
            if type(likelihood_contributions[obs]["uncertainty"]) == float:#usual case, symmetric error
                product_likelihood *= utils.gauss(obsval["value"],likelihood_contributions[obs]["value"],likelihood_contributions[obs]["uncertainty"])
            elif type(likelihood_contributions[obs]["uncertainty"]) == list:#non-symmetric error, use two-sided gaussian
                product_likelihood *= utils.gauss_pm(obsval["value"],likelihood_contributions[obs]["value"],sigma_m = likelihood_contributions[obs]["uncertainty"][0],sigma_p = likelihood_contributions[obs]["uncertainty"][1])

                
    #handle special cases
    #superiso chi2
    chi2 = observables["chi2"]["value"]
    ndf = observables["chi2_ndf"]["value"]
    gamma = math.gamma(float(ndf)/2)
    coeff = pow(chi2,(float(ndf)/2)-1)/((pow(2,float(ndf)/2))*gamma)
    product_likelihood*= (coeff*math.exp(-chi2/2))
    return product_likelihood
def make_decision(candidate_point,prev_l):
    new_point = {}
    L = get_likelihood(candidate_point)
    if L == 0:
        return 0
    R = L/prev_l
    u = random.uniform(0,1)
    if R>u:
        return L
    else:
        return -1
    #make an MCMC decision
