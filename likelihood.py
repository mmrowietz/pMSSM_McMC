import utils
import random
import math

likelihood_contributions ={}

#https://pdg.lbl.gov/2019/tables/rpp2019-sum-quarks.pdf                                                        
likelihood_contributions["mtop"] = {"value":173.1,"uncertainty":0.9}
likelihood_contributions["mbottom"] = {"value":4.18,"uncertainty":[0.02,0.03]}


likelihood_contributions["alpha_s"] = {"value":0.1181,"uncertainty":0.0011}

likelihood_contributions["mhiggs"] = {"value":125.26}#good practise would be to cite from where the value is taken

#https://pdg.lbl.gov/2019/listings/rpp2019-list-w-boson.pdf
likelihood_contributions["mW"] = {"value":80.379,"uncertainty":0.012}

#superiso gives only the SUSY contribution to amu, thus it must be compared against delta amu

#???
likelihood_contributions["a_muon"]={"value":26.8E-10,"uncertainty":4.3E-10}


# https://pdg.lbl.gov/2020/listings/rpp2020-list-B-plus-minus.pdf
likelihood_contributions["BR_B_to_tau_nu"]={"value":1.09E-04,"uncertainty":0.24E-04}

#https://pdg.lbl.gov/2020/listings/rpp2020-list-D-star-2010-plus-minus.pdf
likelihood_contributions["BR_Ds_to_tau_nu"]={"value":5.48E-02,"uncertainty":0.23E-02}
likelihood_contributions["BR_Ds_to_mu_nu"]={"value":5.49E-03,"uncertainty":0.16E-03}



#likelihood_contributions["drho"]={"value":3.9E-04,"uncertainty":1.9E-04}#http://pdg.lbl.gov/2018/reviews/rpp2018-rev-standard-model.pdf


#likelihood_contributions["Delta0_B_to_K_gamma"] = {}
#https://pdg.lbl.gov/2020/listings/rpp2020-list-B-zero.pdf
likelihood_contributions["BR_B0_K0star_gamma"] = {"value":1.24E-05,"uncertainty":0.24E-05}

#https://pdg.lbl.gov/2020/listings/rpp2020-list-Bs-zero.pdf
likelihood_contributions["BR_Bs_to_mu_mu"] = {"value":3.0E-09,"uncertainty":0.4E-09}


#likelihood_contributions["BR_Bd_to_mu_mu"] = {"value":1.1E-10,"uncertainty":[1.3E-10,1.4E-10]}

#likelihood_contributions["BR_b_to_s_mu_mu"] = {}
#likelihood_contributions["BR_b_to_s_e_e"] = {}
#likelihood_contributions["BR_b_to_s_gamma"] = {}

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

    # higgs signals chi2


    # higgs bounds chi2
    
    return product_likelihood


def make_decision(candidate_point,prev_l):
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
