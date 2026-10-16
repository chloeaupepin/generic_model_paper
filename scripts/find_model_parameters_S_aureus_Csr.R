library(dplyr)
library(here)
library(deSolve)

source(here::here("scripts", "model.R"))

simulate_eq_without_vaccine <- function(params) {
  out <- ode(
    y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
          Sv = 0,Csv = 0,Csrv =0,Crv = 0,
          Isnv = 0,Isrnv = 0,Irnv = 0,
          Isv = 0,Isrv = 0,Irv = 0,
          inccumIsnv=0,inccumIsrnv=0, inccumIrnv=0, inccumIsv=0,inccumIsrv=0, inccumIrv=0, inccumSelectionOfResistantForNonVaccinated=0, inccumSelectionOfResistantForVaccinated=0),
    times = seq(from=0,to=1000000,by=100),
    func = SCISsrV.model_3,
    parms = c(params,list(
      vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
      vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0
    ))  
  ) %>% as.data.frame()
  
  out
}


# Compute error between equilibrium value and target
eval_model <- function(x) {
  
  names(x) <- params_to_find
  
  all_params <- c(x,params_fixed_values, list("ar" = x[["as"]]))
  
  out <- simulate_eq_without_vaccine(all_params)
  last <- out[nrow(out), ]
  
  incidenceI <- (last$incidenceIsnv + last$incidenceIsv) +
    (last$incidenceIrnv + last$incidenceIrv)
  
  C  <- last$Csnv +last$Csrnv + last$Crnv
  
  Cr <- last$Crnv
  
  (incidenceI - valeur_cible_I)^2 +
    (C / 100000 - valeur_cible_C)^2 +
    (Cr / C - valeur_cible_Cr)^2
}

check_result <- function(x) {
  
  names(x) <- params_to_find
  
  all_params <- c(x,params_fixed_values, list("ar" = x[["as"]]))
  
  out <- simulate_eq_without_vaccine(all_params)
  last <- out[nrow(out), ]
  
  incidenceI <- (last$incidenceIsnv + last$incidenceIsv) +
    (last$incidenceIrnv + last$incidenceIrv)
  
  C  <- last$Csnv +last$Csrnv + last$Crnv 
  
  Cr <- last$Crnv
  
  
  return(c(incidenceI, C/100000, Cr/C ) )
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
                 "inf_dur_with_abx_r" = 14)

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
                           prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
                           time_until_decolo_by_bystander_ATB = 10,
                           prob_minority_strain_when_colonised = 0.1,
                           eps = 1,
                           prob_specific_exposure = S_aureus$prop_inf_abx_treated_s,
                           time_until_decolo_by_specific_ATB = S_aureus$inf_dur_with_abx_s,
                           prob_minority_strain_when_infected = 0.1,
                           prob_specific_exposure_r = S_aureus$prop_inf_abx_treated_r,
                           time_until_decolo_by_specific_ATB_r = S_aureus$inf_dur_with_abx_r,
                           thetasr = 0,
                           thetars = 0,
                           sigmaC1 = 1,
                           sigmaI1 = 1
)

valeur_cible_I = S_aureus$inc_inf
valeur_cible_C = S_aureus$prop_col
valeur_cible_Cr = S_aureus$prop_col_res_carriers

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
lower <- c(betaC = 0.01, f = 0.92, as = 0.000001)
upper <- c(betaC = 0.02,  f = 0.96,   as = 0.00001)

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
    objective = eval_model,
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
check <- check_result(best$par)
cat("Error =",eval_model(best$par))

#### Save results #### 

S_aureus_params = as.list(c(best$par, params_fixed_values, list("ar" = best$par[["as"]])))

write.csv(S_aureus_params,file = here("files","S_aureus_params_eu_with_J01CR.csv"),row.names=FALSE)
