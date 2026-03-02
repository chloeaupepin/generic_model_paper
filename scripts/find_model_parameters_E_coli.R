
library(dplyr)
library(here)


# Compute error between equilibrium value and target
eval_model_coli <- function(x) {
  # transmission by infected individuals is allowed, antibiotic flux are computed before
  names(x) <- params_to_find
  with(as.list(c(x, params_fixed_values,list("ar" = x[["as"]]), list("betaI" = x[["betaC"]]))), { 
    #print(all_params)
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
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
check_result_coli <- function(x) {
  names(x) <- params_to_find
  with(as.list(c(x, params_fixed_values,list("ar" = x[["as"]]), list("betaI" = x[["betaC"]]))), { 
    #print(all_params)
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
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

params_to_find = c("betaC","f","as")

params_fixed_values = list(time_until_recovery_without_ATB_s = 28,
                           time_until_recovery_without_ATB_r = 28,
                           dps = 365,
                           dpr = 365,
                           prob_bystander_exposure = (0.96+0.77)/1000,
                           time_until_decolo_by_bystander_ATB = 10,
                           prob_minority_strain_when_colonised = 0.1,
                           eps = 1,
                           prob_specific_exposure = 0.8,
                           prob_specific_exposure_r = 0.8,
                           # gammaA_S = prop_first_line/temps_first_line + prop_second_line/temps_second_line*(1-prop_min) 
                           gammaA_s = 0.75/2+0.05/7*(1-0.1), 
                           # ANTIBIO SECOND LINE DOIT SELECTIONNER POUR LA RESISTANCE
                           # psiA = prop_second_line/temps_second_line*prop_min 
                           psiA = 0.05/7*0.1, 
                           # gammaA_R = prop_first_line/temps_first_line + prop_second_line/temps_second_line
                           gammaA_r = 0.75/2+0.05/9, 
                           prob_minority_strain_when_infected = 0.1,
                           thetasr = 0,
                           thetars = 0
)

valeur_cible_I = 16
valeur_cible_C = 0.95
valeur_cible_Cr = 0.06

#### Optimise parameters ####

params_to_find_init_values = list(betaC = 0.1, f = 0.87, as = 0.000004)
#eval_model_coli(params_to_find_init_values)

res <- nlminb(start = params_to_find_init_values,
              objective = eval_model_coli,
              lower = c(0.01,0.85,0.000001),#c(0.1,0.95,0.000001,0.000001), #
              upper = c(0.5,0.99,0.0005))#c(0.5,1,0.0005,0.0005)) #
print(res$par)
check <- check_result_coli(res$par)
eval_model_coli(res$par)

#### Save results #### 

E_coli_params = as.list(c(res$par, params_fixed_values, list("ar" = res$par[["as"]]), list("betaI" = res$par[["betaC"]])))

write.csv(E_coli_params,file = here::here("case_studies","data","E_coli_params_primavera4.csv"),row.names=FALSE)

