library(dplyr)
library(here)


# Compute error between equilibrium value and target

eval_model_coli <- function(x) {
  # transmission by infected individuals is allowed, antibiotic flux are computed before
  names(x) <- params_to_find
  with(as.list(c(x, params_fixed_values,list("ar" = x[["as"]]), list("betaI" = x[["betaC"]]))), { 
    #print(all_params)
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    # gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB
    
    # gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    
    etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    # Sigma 2 
    sigmaC2 = lambdaA
    sigmaI2 = lambdaA_I + gammaA_s
    
    His = etas +lambdaA_I + gammaA_s + sigmaI1
    Hisr = etas + sigmaI2
    
    N = 100000
    X = (sigmaC1 + etas*sigmaI1/Hisr*as/His)/(1/dps + as + sigmaC2 - etas*as/Hisr)
    
    eq_Snv = N/(betaC*(1+X) + betaI*(as/His+ as/Hisr*X + sigmaI1/Hisr*as/His))*(1/dps + lambdaA + sigmaC1 + as - etas*as/His)
    
    E1 = f*eq_Snv/N*(betaC + betaI*ar/(etar+gammaA_r))-1/dpr - ar + etar*ar/(etar+gammaA_r)
    E2 = (betaI*f*eq_Snv/N +etar)* sigmaI2/(etar+gammaA_r) * (as/Hisr*X + sigmaI1/Hisr*as/His)+sigmaC2*X
    
    E = E1/E2
    
    G = 1+ ar/(etar+gammaA_r)-E*(1+X+as/His+(as/Hisr*X + sigmaI1/Hisr*as/His)*(1+sigmaI2/(etar+gammaA_r)))
    
    eq_Cr = (N-eq_Snv)/G
    eq_Cs = -E*eq_Cr
    eq_Csr = X*eq_Cs
    
    
    error = (as*eq_Cs + as*eq_Csr + ar*eq_Cr - valeur_cible_I)^2+
      ((eq_Cs + eq_Csr + eq_Cr)/100000 - valeur_cible_C)^2 + ((eq_Cr/(eq_Cs + eq_Cr+ eq_Csr) - valeur_cible_Cr))^2
    
    # Erreur au carré
    return(error)
  })
}


# Compute equilibrium value (possibility to compare afterwarsds with target value)

check_result_coli <- function(x) {
  names(x) <- params_to_find
  with(as.list(c(x, params_fixed_values,list("ar" = x[["as"]]), list("betaI" = x[["betaC"]]))), { 
    #print(all_params)
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    # gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB
    # 
    # gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    
    etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    # Sigma 2 
    sigmaC2 = lambdaA
    sigmaI2 = lambdaA_I + gammaA_s
    
    His = etas +lambdaA_I + gammaA_s + sigmaI1
    Hisr = etas + sigmaI2
    
    N = 100000
    X = (sigmaC1 + etas*sigmaI1/Hisr*as/His)/(1/dps + as + sigmaC2 - etas*as/Hisr)
    
    eq_Snv = N/(betaC*(1+X) + betaI*(as/His+ as/Hisr*X + sigmaI1/Hisr*as/His))*(1/dps + lambdaA + sigmaC1 + as - etas*as/His)
    
    E1 = f*eq_Snv/N*(betaC + betaI*ar/(etar+gammaA_r))-1/dpr - ar + etar*ar/(etar+gammaA_r)
    E2 = (betaI*f*eq_Snv/N +etar)* sigmaI2/(etar+gammaA_r) * (as/Hisr*X + sigmaI1/Hisr*as/His)+sigmaC2*X
    
    E = E1/E2
    
    G = 1+ ar/(etar+gammaA_r)-E*(1+X+as/His+(as/Hisr*X + sigmaI1/Hisr*as/His)*(1+sigmaI2/(etar+gammaA_r)))
    
    eq_Cr = (N-eq_Snv)/G
    eq_Cs = -E*eq_Cr
    eq_Csr = X*eq_Cs
    
    
    cat("=== Target check ===\n")
    
    cat(sprintf("Target 1: %.4f (expected: %.4f) \n",
                as*eq_Cs + as*eq_Csr + ar*eq_Cr, valeur_cible_I))
    
    cat(sprintf("Target 2: %.4f (expected: %.4f)\n",
                (eq_Cs + eq_Cr + eq_Csr)/100000, valeur_cible_C))
    
    cat(sprintf("Target 3: %.4f (expected: %.4f)\n",
                eq_Cr/(eq_Cs + eq_Cr + eq_Csr), valeur_cible_Cr))
    
    result = c(as*eq_Cs + as*eq_Csr + ar*eq_Cr, (eq_Cs + eq_Cr + eq_Csr)/100000, eq_Cr/(eq_Cs + eq_Cr + eq_Csr))
    
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
                           gammaA_s = 0.6/2+0.2/7*(1-0.1), 
                           # ANTIBIO SECOND LINE DOIT SELECTIONNER POUR LA RESISTANCE
                           # psiA = prop_second_line/temps_second_line*prop_min 
                           psiA = 0.2/7*0.1, 
                           # gammaA_R = prop_first_line/temps_first_line + prop_second_line/temps_second_line
                           gammaA_r = 0.6/2+0.2/9, 
                           prob_minority_strain_when_infected = 0.1,
                           thetasr = 0,
                           thetars = 0,
                           sigmaC1 = 1e-4, 
                           sigmaI1 = 1e-4
)

valeur_cible_I = 16
valeur_cible_C = 0.95
valeur_cible_Cr = 0.06

#### Optimise parameters ####

# params_to_find_init_values = list(betaC = 0.015, f = 0.85, as = 0.000004)
# eval_model(params_to_find_init_values)
# 
# res <- nlminb(start = params_to_find_init_values,
#               objective = eval_model,
#               lower = c(0.005,0.8,0.000001),
#               upper = c(0.05,1,0.00001))
# print(res$par)
# check <- check_result(res$par)
# eval_model(res$par)

#### Optimise parameters with multiple start ####

# Bounds
lower <- c(betaC = 0.01, f = 0.8, as = 0.000001)
upper <- c(betaC = 0.5,  f = 1,   as = 0.0005)

# Number of points by dimension
n_grid <- 10

# Grid creation
grid <- expand.grid(
  betaC = seq(lower["betaC"], upper["betaC"], length.out = n_grid),
  f     = seq(lower["f"],     upper["f"],     length.out = n_grid),
  as    = seq(lower["as"],    upper["as"],    length.out = n_grid)
)

# Total number of combinations
n_starts <- nrow(grid)

# Storage
results <- vector("list", n_starts)

for (i in 1:n_starts) {
  start_vals <- as.numeric(grid[i, ])
  names(start_vals) <- names(lower)
  
  res <- nlminb(
    start = start_vals,
    objective = eval_model_coli,
    lower = lower,
    upper = upper
  )
  
  results[[i]] <- list(
    par = res$par,
    obj = res$objective,
    convergence = res$convergence,
    start = start_vals
  )
}

results_df <- data.frame(do.call(rbind, lapply(results, function(x) {
  c(x$start, x$par, obj = x$obj, conv = x$convergence)
})))

# Best parameters
best_index <- which.min(sapply(results, function(x) x$obj))
best <- results[[best_index]]

cat("Optimised parameters:\n", paste0(names(best$par)," = ", best$par, "\n"))

# Check
check <- check_result_coli(best$par)
cat("Error =",eval_model_coli(best$par))

#### Save results #### 

E_coli_params = as.list(c(best$par, params_fixed_values, list("ar" = best$par[["as"]]), list("betaI" = best$par[["betaC"]])))

write.csv(E_coli_params,file = here::here("files","E_coli_params_2.csv"),row.names=FALSE)
