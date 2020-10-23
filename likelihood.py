import utils
import random



likelihood_contributions ={}
likelihood_contributions["mtop"] = {"value":173.1,"uncertainty":0.9}
likelihood_contributions["mbottom"] = {"value":4.18,"uncertainty":[0.03,0.04]}
likelihood_contributions["alpha_s"] = {"value":0.1181,"uncertainty":0.0011}
likelihood_contributions["mhiggs"] = {"value":125.26}
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
            exit()
        if "uncertainty" in obsval.keys(): #higgs-type case: center gauss on theory value with width= theory uncertainty, evaluate at experimental value
            product_likelihood *= utils.gauss(likelihood_contributions[obs]["value"],obsval["value"],0.5*obsval["uncertainty"])
        else:
            if type(likelihood_contributions[obs]["uncertainty"]) == float:#usual case, symmetric error
                product_likelihood *= utils.gauss(obsval["value"],likelihood_contributions[obs]["value"],likelihood_contributions[obs]["uncertainty"])
            elif type(likelihood_contributions[obs]["uncertainty"]) == list:#non-symmetric error, use two-sided gaussian
                product_likelihood *= utils.gauss_pm(obsval["value"],likelihood_contributions[obs]["value"],sigma_m = likelihood_contributions[obs]["uncertainty"][0],sigma_p = likelihood_contributions[obs]["uncertainty"][1])


    #handle special cases
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
