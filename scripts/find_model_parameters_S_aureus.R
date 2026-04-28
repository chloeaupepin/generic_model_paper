library(dplyr)
library(here)


# Compute error between equilibrium value and target

eval_model <- function(x) {
  names(x) <- params_to_find
  with(as.list(c(x, params_fixed_values,list("ar" = x[["as"]]))), { 
    #print(all_params)
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
    gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*(1-prob_minority_strain_when_infected)
    psiA = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*prob_minority_strain_when_infected
    
    gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    etas = (1-prob_specific_exposure) * 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    N = 100000
    His = etas + lambdaA_I +  gammaA_s+ psiA + phiA_I
    eq_Snv = N/(betaC + betaI*as/His)*(1/dps + lambdaA + as - etas*as/His + phiA)
    E1 = f*eq_Snv/N*(betaC + betaI*ar/(etar+gammaA_r))-1/dpr - ar + etar*ar/(etar+gammaA_r)
    E2 = (betaI*f*eq_Snv/N +etar)*(psiA + phiA_I)/(etar+gammaA_r)*as/His +phiA
    E = E1/E2
    G = 1+ ar/(etar+gammaA_r)-E*(1+as/His*(1+(psiA + phiA_I)/(etar+gammaA_r)))
    eq_Cr = (N-eq_Snv)/G
    eq_Cs = -E*eq_Cr
    
    
    error = (as*eq_Cs + ar*eq_Cr - valeur_cible_I)^2+
      ((eq_Cs + eq_Cr)/100000 - valeur_cible_C)^2 + ((eq_Cr/(eq_Cs + eq_Cr) - valeur_cible_Cr))^2
    
    # Erreur au carré
    return(error)
  })
}


# Compute equilibrium value (possibility to compare afterwarsds with target value)

check_result <- function(x) {
  names(x) <- params_to_find
  with(as.list(c(x, params_fixed_values,list("ar" = x[["as"]]))), { 
    #print(all_params)
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
    gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*(1-prob_minority_strain_when_infected)
    psiA = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*prob_minority_strain_when_infected
    
    gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    etas = (1-prob_specific_exposure) * 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    N = 100000
    His = etas + lambdaA_I +  gammaA_s+ psiA + phiA_I
    eq_Snv = N/(betaC + betaI*as/His)*(1/dps + lambdaA + as - etas*as/His + phiA)
    E1 = f*eq_Snv/N*(betaC + betaI*ar/(etar+gammaA_r))-1/dpr - ar + etar*ar/(etar+gammaA_r)
    E2 = (betaI*f*eq_Snv/N +etar)*(psiA + phiA_I)/(etar+gammaA_r)*as/His +phiA
    E = E1/E2
    G = 1+ ar/(etar+gammaA_r)-E*(1+as/His*(1+(psiA + phiA_I)/(etar+gammaA_r)))
    eq_Cr = (N-eq_Snv)/G
    eq_Cs = -E*eq_Cr
    
    cat("=== Target check ===\n")
    
    cat(sprintf("Target 1: %.4f (expected: %.4f) \n",
                as*eq_Cs + ar*eq_Cr, valeur_cible_I))
    
    cat(sprintf("Target 2: %.4f (expected: %.4f)\n",
                (eq_Cs + eq_Cr)/100000, valeur_cible_C))
    
    cat(sprintf("Target 3: %.4f (expected: %.4f)\n",
                eq_Cr/(eq_Cs + eq_Cr), valeur_cible_Cr))
    
    result = c(as*eq_Cs + ar*eq_Cr, (eq_Cs + eq_Cr)/100000, eq_Cr/(eq_Cs + eq_Cr))
    
    return(result)
  })
}

#### Define parameter values ####

S_aureus <- list("prop_col" = 0.3,
                 "prop_col_res_carriers" = 0.1,
                 "inc_inf" = 37.9/365,
                 #"inc_inf_s" = 22.7/365, #5,
                 #"inc_inf_r" = 0.18/365, #5,
                 "carriage_dur_s" = 98,
                 "carriage_dur_r" = 98,
                 "inf_dur_without_abx_s" = 50,
                 "inf_dur_without_abx_r" = 50,
                 "prop_inf_abx_treated_s" = 1, #O.75
                 "prop_inf_abx_treated_r" = 1, #O.75
                 "inf_dur_with_abx_s" = 14,
                 "inf_dur_with_abx_r" = 28)

params_to_find = c("betaC","f","as")

params_fixed = c("betaI","etas","etar","dps","dpr","dir",
                 "prob_bystander_exposure", "time_until_decolo_by_bystander_ATB",
                 "prob_minority_strain_when_colonised",
                 "eps","prob_specific_exposure",
                 "time_until_decolo_by_specific_ATB",
                 "prob_minority_strain_when_infected",
                 "thetasr",
                 "thetars"
)

params_fixed_values = list(betaI = 0, 
                           time_until_recovery_without_ATB_s = S_aureus$inf_dur_without_abx_s,
                           time_until_recovery_without_ATB_r = S_aureus$inf_dur_without_abx_r,
                           dps = S_aureus$carriage_dur_s,
                           dpr = S_aureus$carriage_dur_r,
                           prob_bystander_exposure = 0.06/1000,
                           time_until_decolo_by_bystander_ATB = 10,
                           prob_minority_strain_when_colonised = 0.1,
                           eps = 1,
                           prob_specific_exposure = S_aureus$prop_inf_abx_treated_s,
                           time_until_decolo_by_specific_ATB = S_aureus$inf_dur_with_abx_s,
                           prob_minority_strain_when_infected = 0.1,
                           prob_specific_exposure_r = S_aureus$prop_inf_abx_treated_r,
                           time_until_decolo_by_specific_ATB_r = S_aureus$inf_dur_with_abx_r,
                           thetasr = 0,
                           thetars = 0
)

valeur_cible_I = S_aureus$inc_inf
valeur_cible_C = S_aureus$prop_col
valeur_cible_Cr = S_aureus$prop_col_res_carriers

#### Optimise parameters ####

params_to_find_init_values = list(betaC = 0.015, f = 0.99, as = 0.000004)
#eval_model(params_to_find_init_values)

res <- nlminb(start = params_to_find_init_values,
              objective = eval_model,
              lower = c(0.005,0.95,0.000001),
              upper = c(0.05,1,0.00001))
print(res$par)
check <- check_result(res$par)
eval_model(res$par)

#### Save results #### 

S_aureus_params = as.list(c(res$par, params_fixed_values, list("ar" = res$par[["as"]])))

write.csv(S_aureus_params,file = here("files","S_aureus_params_eu.csv"),row.names=FALSE)
