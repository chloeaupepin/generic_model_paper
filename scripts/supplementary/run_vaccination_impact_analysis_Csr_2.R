library(purrr)
library(tidyr)
library(dplyr)
library(deSolve)
library(furrr)
library(stringr)

source(here::here("scripts", "model.R"))
source(here::here("scripts","utils", "functions.R"))

compute_equilibrium <- function(data,population_size = 100000){
  
  data %>%
    mutate(
      lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB,
      
      lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB,
      
      # gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB,
      # 
      # gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r,
      
      etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s,
      etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r,
      
      # Sigma 2 
      sigmaC2 = lambdaA,
      sigmaI2 = lambdaA_I + gammaA_s,
      
      His = etas +lambdaA_I + gammaA_s + sigmaI1,
      Hisr = etas + sigmaI2,
      
      N = population_size,
      X = (sigmaC1 + etas*sigmaI1/Hisr*as/His)/(1/dps + as + sigmaC2 - etas*as/Hisr),
      
      eq_Snv = N/(betaC*(1+X) + betaI*(as/His+ as/Hisr*X + sigmaI1/Hisr*as/His))*(1/dps + lambdaA + sigmaC1 + as - etas*as/His),
      
      E1 = f*eq_Snv/N*(betaC + betaI*ar/(etar+gammaA_r))-1/dpr - ar + etar*ar/(etar+gammaA_r),
      E2 = (betaI*f*eq_Snv/N +etar)* sigmaI2/(etar+gammaA_r) * (as/Hisr*X + sigmaI1/Hisr*as/His)+sigmaC2*X,
      
      E = E1/E2,
      
      G = 1+ ar/(etar+gammaA_r)-E*(1+X+as/His+(as/Hisr*X + sigmaI1/Hisr*as/His)*(1+sigmaI2/(etar+gammaA_r))),
      
      eq_Crnv = (N-eq_Snv)/G,
      eq_Csnv = -E*eq_Crnv,
      eq_Csrnv = X*eq_Csnv,
      eq_Isnv = as/His*eq_Csnv,
      eq_Isrnv = as/Hisr*eq_Csrnv + sigmaI1/Hisr*eq_Isnv,
      eq_Irnv = ar/(etar+gammaA_r)*eq_Crnv + sigmaI2/(etar+gammaA_r)*eq_Isrnv
    )
}

# Function to run model one year without vaccine

simulate_1y_without_vaccine <- function(params) {
  
  out <- ode(
    y = 
      c(Snv = params$eq_Snv,Csnv = params$eq_Csnv,Csrnv = params$eq_Csrnv,Crnv = params$eq_Crnv,
        Sv = 0,Csv = 0, Csrv = 0, Crv =0, 
        Isnv = params$eq_Isnv ,Isrnv = params$eq_Isrnv ,Irnv = params$eq_Irnv , Isv = 0,Isrv = 0, Irv =0,
        inccumIsnv=0,inccumIsrnv=0, inccumIrnv=0, inccumIsv=0,inccumIsrv=0, inccumIrv=0, inccumSelectionOfResistantForNonVaccinated=0, inccumSelectionOfResistantForVaccinated=0),
    times = seq(from=0,to=370,by=10),
    func = SCISsrV.model_2,
    parms = c(params, vaccine_params_zero())  
  ) %>% as.data.frame()
  
  out %>%
    slice_tail(n=1) %>% 
    select(Snv, Csnv, Csrnv, Crnv, Isnv,Isrnv, Irnv,
           inccumIsnv,inccumIsrnv, inccumIrnv, inccumIsv,inccumIsrv, inccumIrv,
           inccumSelectionOfResistantForNonVaccinated,inccumSelectionOfResistantForVaccinated)
}

# Function to run model one year with vaccine

simulate_1y_with_vaccine <- function(params) {
  out <- ode(
    y = 
      c(Snv = params$eq_Snv * (1-params$Vperc),Csnv = params$eq_Csnv * (1-params$Vperc),Csrnv = params$eq_Csrnv * (1-params$Vperc),Crnv = params$eq_Crnv * (1-params$Vperc),
        Sv = params$eq_Snv * (params$Vperc),Csv = params$eq_Csnv * (params$Vperc),Csrv = params$eq_Csrnv * (params$Vperc), Crv =params$eq_Crnv * (params$Vperc), 
        Isnv = params$eq_Isnv * (1-params$Vperc),Isrnv = params$eq_Isrnv * (1-params$Vperc),Irnv = params$eq_Irnv * (1-params$Vperc), Isv = params$eq_Isnv * (params$Vperc),Isrv = params$eq_Isrnv * (params$Vperc), Irv =params$eq_Irnv * (params$Vperc),
        inccumIsnv=0,inccumIsrnv=0, inccumIrnv=0, inccumIsv=0,inccumIsrv=0, inccumIrv=0, inccumSelectionOfResistantForNonVaccinated=0, inccumSelectionOfResistantForVaccinated=0),
    times = seq(from=0,to=370,by=10),
    func = SCISsrV.model_2,
    parms = c(params)  
  ) %>% as.data.frame()
  
  out %>%
    slice_tail(n=1) %>% 
    select(Snv, Csnv,Csrnv, Crnv, Isnv,Isrnv, Irnv, Sv, Csv,Csrv, Crv, Isv,Isrv, Irv,
           inccumIsnv, inccumIsrnv, inccumIrnv, inccumIsv,inccumIsrv, inccumIrv,
           inccumSelectionOfResistantForNonVaccinated,inccumSelectionOfResistantForVaccinated)
}

compute_outputs <- function(sim, population_size = 100000, factor = 100000){
  # Add cumulative incidence of selection to cumulative incidence of resistant infections
  results <- sim %>%
    mutate(res_1y_wov_inccumIrnv_sel = res_1y_wov_inccumIrnv + res_1y_wov_inccumSelectionOfResistantForNonVaccinated,
           res_1y_wv_inccumIrnv_sel = res_1y_wv_inccumIrnv + res_1y_wv_inccumSelectionOfResistantForNonVaccinated,
           res_1y_wv_inccumIrv_sel = res_1y_wv_inccumIrv + res_1y_wv_inccumSelectionOfResistantForVaccinated)
  
  
  # Add cumulative incidence of total infections
  results <- results %>%
    mutate(res_1y_wov_inccumI = res_1y_wov_inccumIrnv_sel+res_1y_wov_inccumIsnv+res_1y_wov_inccumIsrnv,
           res_1y_wv_inccumI_non_vaccinated = res_1y_wv_inccumIrnv_sel + res_1y_wv_inccumIsnv+ res_1y_wv_inccumIsrnv,
           res_1y_wv_inccumI_vaccinated = res_1y_wv_inccumIrv_sel + res_1y_wv_inccumIsv+ res_1y_wv_inccumIsrv,
           res_1y_wv_inccumI = res_1y_wv_inccumI_non_vaccinated + res_1y_wv_inccumI_vaccinated)%>%
    mutate(res_1y_wov_inccumI_old = res_1y_wov_inccumIrnv+res_1y_wov_inccumIsnv+res_1y_wov_inccumIsrnv,
           res_1y_wv_inccumI_old_non_vaccinated = res_1y_wv_inccumIrnv + res_1y_wv_inccumIsnv+ res_1y_wv_inccumIsrnv,
           res_1y_wv_inccumI_old_vaccinated = res_1y_wv_inccumIrv + res_1y_wv_inccumIsv+ res_1y_wv_inccumIsrv,
           res_1y_wv_inccumI_old = res_1y_wv_inccumI_old_non_vaccinated + res_1y_wv_inccumI_old_vaccinated)
  
  # Compute reduction change in cumulative incidences
  results<- results %>%
    mutate(prc_red_inccumIs = prc_red(res_1y_wov_inccumIsnv + res_1y_wov_inccumIsrnv,res_1y_wv_inccumIsnv + res_1y_wv_inccumIsv + res_1y_wv_inccumIsrnv + res_1y_wv_inccumIsrv),
           prc_red_inccumIs_non_vaccinated = prc_red((res_1y_wov_inccumIsnv + res_1y_wov_inccumIsrnv)*(1-Vperc),res_1y_wv_inccumIsnv + res_1y_wv_inccumIsrnv ),
           prc_red_inccumIs_vaccinated = prc_red((res_1y_wov_inccumIsnv + res_1y_wov_inccumIsrnv)*Vperc, res_1y_wv_inccumIsv + res_1y_wv_inccumIsrv),
           prc_red_inccumIr = prc_red(res_1y_wov_inccumIrnv_sel,res_1y_wv_inccumIrnv_sel + res_1y_wv_inccumIrv_sel),
           prc_red_inccumIr_non_vaccinated = prc_red(res_1y_wov_inccumIrnv_sel*(1-Vperc),res_1y_wv_inccumIrnv_sel ),
           prc_red_inccumIr_vaccinated = prc_red(res_1y_wov_inccumIrnv_sel*Vperc, res_1y_wv_inccumIrv_sel),
           prc_red_inccumI = prc_red(res_1y_wov_inccumI, res_1y_wv_inccumI),
           prc_red_inccumI_non_vaccinated = prc_red(res_1y_wov_inccumI*(1-Vperc), res_1y_wv_inccumI_non_vaccinated),
           prc_red_inccumI_vaccinated = prc_red(res_1y_wov_inccumI*Vperc, res_1y_wv_inccumI_vaccinated)
    )
  
  # Compute proportion of resistant infections among infections in cumulative incidences
  results <- results %>%
    mutate(res_1y_wov_prop_inccumIr_sel = (res_1y_wov_inccumIrnv_sel)/res_1y_wov_inccumI,
           res_1y_wv_prop_inccumIr_sel = (res_1y_wv_inccumIrnv_sel + res_1y_wv_inccumIrv_sel)/res_1y_wv_inccumI,
           res_1y_wv_prop_inccumIr_sel_non_vaccinated = (res_1y_wv_inccumIrnv_sel)/res_1y_wv_inccumI_non_vaccinated,
           res_1y_wv_prop_inccumIr_sel_vaccinated = (res_1y_wv_inccumIrv_sel)/res_1y_wv_inccumI_vaccinated
    )%>%
    mutate(res_1y_wov_prop_inccumIr = (res_1y_wov_inccumIrnv)/(res_1y_wov_inccumIrnv+res_1y_wov_inccumIsnv),
           res_1y_wv_prop_inccumIr = (res_1y_wv_inccumIrnv + res_1y_wv_inccumIrv)/(res_1y_wv_inccumIrnv + res_1y_wv_inccumIsnv + res_1y_wv_inccumIrv + res_1y_wv_inccumIsv),
           res_1y_wv_prop_inccumIr_non_vaccinated = (res_1y_wv_inccumIrnv)/( res_1y_wv_inccumIrnv + res_1y_wv_inccumIsnv),
           res_1y_wv_prop_inccumIr_vaccinated = (res_1y_wv_inccumIrv)/(res_1y_wv_inccumIrv + res_1y_wv_inccumIsv)
    )
  
  
  # Compute proportion of resistant colonisation among colonisation in prevalences at the end of the year
  results <- results %>%
    mutate(res_1y_wv_prevC =( res_1y_wv_Crnv + res_1y_wv_Csnv +res_1y_wv_Csrnv + res_1y_wv_Crv  + res_1y_wv_Csv+ res_1y_wv_Csrv)/population_size) %>%
    mutate(res_1y_wov_prop_prevCIr = (res_1y_wov_Crnv+res_1y_wov_Irnv)/(res_1y_wov_Crnv + res_1y_wov_Irnv + res_1y_wov_Csnv + res_1y_wov_Isnv + res_1y_wov_Csrnv + res_1y_wov_Isrnv),
           res_1y_wv_prop_prevCIr = (res_1y_wv_Crnv+res_1y_wv_Irnv + res_1y_wv_Crv+res_1y_wv_Irv)/(res_1y_wv_Crnv + res_1y_wv_Irnv + res_1y_wv_Csnv + res_1y_wv_Isnv + res_1y_wv_Csrnv + res_1y_wv_Isrnv+ res_1y_wv_Crv + res_1y_wv_Irv + res_1y_wv_Csv + res_1y_wv_Isv + res_1y_wv_Csrv + res_1y_wv_Isrv),
           res_1y_wv_prop_prevCIr_non_vaccinated = (res_1y_wv_Crnv+res_1y_wv_Irnv)/(res_1y_wv_Crnv + res_1y_wv_Irnv + res_1y_wv_Csnv + res_1y_wv_Isnv+ res_1y_wv_Csrnv + res_1y_wv_Isrnv),
           res_1y_wv_prop_prevCIr_vaccinated = (res_1y_wv_Crv+res_1y_wv_Irv)/(res_1y_wv_Crv + res_1y_wv_Irv + res_1y_wv_Csv + res_1y_wv_Isv  + res_1y_wv_Csrv + res_1y_wv_Isrv),
           
           res_1y_wov_prop_prevCr = (res_1y_wov_Crnv)/(res_1y_wov_Crnv + res_1y_wov_Csnv + res_1y_wov_Csrnv),
           res_1y_wv_prop_prevCr = (res_1y_wv_Crnv + res_1y_wv_Crv)/(res_1y_wv_Crnv + res_1y_wv_Csnv + res_1y_wv_Csrnv + res_1y_wv_Crv  + res_1y_wv_Csv  + res_1y_wv_Csrv ),
           res_1y_wv_prop_prevCr_non_vaccinated = (res_1y_wv_Crnv)/(res_1y_wv_Crnv + res_1y_wv_Csnv+ res_1y_wv_Csrnv ),
           res_1y_wv_prop_prevCr_vaccinated = (res_1y_wv_Crv)/(res_1y_wv_Crv + res_1y_wv_Csv  + res_1y_wv_Csrv )
    )
  
  # Compute reduction change in proportions
  results<- results %>%
    mutate(prc_red_prop_prevCIr = prc_red(res_1y_wov_prop_prevCIr,res_1y_wv_prop_prevCIr ),
           prc_red_prop_prevCIr_non_vaccinated = prc_red(res_1y_wov_prop_prevCIr,res_1y_wv_prop_prevCIr_non_vaccinated ),
           prc_red_prop_prevCIr_vaccinated = prc_red(res_1y_wov_prop_prevCIr,res_1y_wv_prop_prevCIr_vaccinated ))%>%
    mutate(prc_red_prop_prevCr = prc_red(res_1y_wov_prop_prevCr,res_1y_wv_prop_prevCr ),
           prc_red_prop_prevCr_non_vaccinated = prc_red(res_1y_wov_prop_prevCr,res_1y_wv_prop_prevCr_non_vaccinated ),
           prc_red_prop_prevCr_vaccinated = prc_red(res_1y_wov_prop_prevCr,res_1y_wv_prop_prevCr_vaccinated ))%>%
    mutate(prc_red_prop_inccumIr_sel = prc_red(res_1y_wov_prop_inccumIr_sel,res_1y_wv_prop_inccumIr_sel ),
           prc_red_prop_inccumIr_sel_non_vaccinated = prc_red(res_1y_wov_prop_inccumIr_sel,res_1y_wv_prop_inccumIr_sel_non_vaccinated ),
           prc_red_prop_inccumIr_sel_vaccinated = prc_red(res_1y_wov_prop_inccumIr_sel,res_1y_wv_prop_inccumIr_sel_vaccinated ))%>%
    mutate(prc_red_prop_inccumIr = prc_red(res_1y_wov_prop_inccumIr,res_1y_wv_prop_inccumIr ),
           prc_red_prop_inccumIr_non_vaccinated = prc_red(res_1y_wov_prop_inccumIr,res_1y_wv_prop_inccumIr_non_vaccinated ),
           prc_red_prop_inccumIr_vaccinated = prc_red(res_1y_wov_prop_inccumIr,res_1y_wv_prop_inccumIr_vaccinated ))
  
  # Compute reduction change in prevalences
  results <- results %>%
    mutate(prc_red_prevCr = prc_red(res_1y_wov_Crnv, res_1y_wv_Crnv + res_1y_wv_Crv),
           prc_red_prevCr_non_vaccinated = prc_red(res_1y_wov_Crnv*(1-Vperc), res_1y_wv_Crnv),
           prc_red_prevCr_vaccinated = prc_red(res_1y_wov_Crnv*Vperc, res_1y_wv_Crv))%>%
    mutate(prc_red_prevC = prc_red(res_1y_wov_Crnv+res_1y_wov_Csnv+res_1y_wov_Csrnv,res_1y_wv_Crnv + res_1y_wv_Crv + res_1y_wv_Csnv + res_1y_wv_Csv+ res_1y_wv_Csrnv + res_1y_wv_Csrv),
           prc_red_prevC_non_vaccinated = prc_red((res_1y_wov_Crnv+res_1y_wov_Csnv+res_1y_wov_Csrnv)*(1-Vperc),res_1y_wv_Crnv +res_1y_wv_Csnv+res_1y_wv_Csrnv ),
           prc_red_prevC_vaccinated = prc_red((res_1y_wov_Crnv+res_1y_wov_Csnv+res_1y_wov_Csrnv)*Vperc, res_1y_wv_Crv + res_1y_wv_Csv+ res_1y_wv_Csrv))
  
  
  # Compute averted infections 
  results <- results %>%
    mutate(averted_inccumIs = res_1y_wov_inccumIsnv + res_1y_wov_inccumIsrnv - res_1y_wv_inccumIsnv - res_1y_wv_inccumIsv- res_1y_wv_inccumIsrnv - res_1y_wv_inccumIsrv,
           averted_inccumIr = res_1y_wov_inccumIrnv_sel - res_1y_wv_inccumIrnv_sel - res_1y_wv_inccumIrv_sel,
           averted_inccumI = averted_inccumIs + averted_inccumIr,
           averted_inccumIs_non_vaccinated = (res_1y_wov_inccumIsnv +res_1y_wov_inccumIsrnv )*(1-Vperc) - res_1y_wv_inccumIsnv - res_1y_wv_inccumIsrnv ,
           averted_inccumIr_non_vaccinated = res_1y_wov_inccumIrnv_sel*(1-Vperc) - res_1y_wv_inccumIrnv_sel,
           averted_inccumI_non_vaccinated = averted_inccumIs_non_vaccinated + averted_inccumIr_non_vaccinated,
           averted_inccumIs_vaccinated = (res_1y_wov_inccumIsnv + res_1y_wov_inccumIsrnv)*Vperc - res_1y_wv_inccumIsv - res_1y_wv_inccumIsrv,
           averted_inccumIr_vaccinated = res_1y_wov_inccumIrnv_sel*Vperc - res_1y_wv_inccumIrv_sel,
           averted_inccumI_vaccinated = averted_inccumIs_vaccinated + averted_inccumIr_vaccinated)%>%
    mutate(across(starts_with("averted"), ~ .x * factor/population_size))
  
  return(results)
  
}



n_cores = parallel::detectCores() - 1
plan(multisession, workers = n_cores)


set.seed(123)  # global seed for reproducibility

#### Define settings ####

## Population size
population_size = 100000

## Choose bacteria 
# Comment or uncomment 4 line codes below and rerun everything for each bacteria 
# S_aureus
bacteria = "S_aureus"
folder_name = "S_aureus"
file_name  = "S_aureus_params_eu_with_J01CR_2.csv"
transmission_by_infected = FALSE

# E_coli
bacteria = "E_coli"
folder_name = "E_coli"
file_name = "E_coli_params_2.csv"#"E_coli_params_primavera5.csv"
transmission_by_infected = TRUE # if transmission is allowed, it equals transmission by colonized individuals

## Choose number of bacteria to generate
n <- 500


#### Define bacteria ####
Bacteria_params = read.csv(here::here("files",file_name)) 

# cols_to_noise include bacteria parameters and antibiotic parameters
# the relative fitness f isn't varied due to too much sensitivity of resistant proportion
if(bacteria == "S_aureus"){
  cols_to_noise <- c(	
    "betaC","as","time_until_recovery_without_ATB_s","dps",
    "prob_bystander_exposure",
    "time_until_decolo_by_bystander_ATB",
    "prob_minority_strain_when_colonised",
    "prob_specific_exposure",
    "time_until_decolo_by_specific_ATB",
    "prob_minority_strain_when_infected",
    "prob_specific_exposure_r",
    "time_until_decolo_by_specific_ATB_r"
  )
} else {
  cols_to_noise <- c(	
    "betaC","as","time_until_recovery_without_ATB_s","dps",
    "prob_bystander_exposure",
    "time_until_decolo_by_bystander_ATB",
    "prob_minority_strain_when_colonised",
    "prob_specific_exposure",
    "prob_minority_strain_when_infected",
    "prob_specific_exposure_r",
    "gammaA_s",
    "psiA",
    "gammaA_r"
  )
}

df_bacteria <- add_variability(Bacteria_params, cols_to_noise, n ,0.1) 

df_bacteria <- set_resistant_and_sensitive_param_equal(df_bacteria, c("dpr", "ar", "time_until_recovery_without_ATB_r"))


if(transmission_by_infected){
  df_bacteria <- df_bacteria %>%
    mutate(betaI = betaC)
} else {
  df_bacteria <- df_bacteria %>%
    mutate(betaI = 0)
}


df_bacteria <- df_bacteria %>%
  mutate(prob_specific_exposure_r = pmin(1,prob_specific_exposure_r),
         prob_specific_exposure = pmin(1,prob_specific_exposure),
         prob_minority_strain_when_colonised = pmin(1,prob_minority_strain_when_colonised),
         prob_minority_strain_when_infected = pmin(1,prob_minority_strain_when_infected),
         prob_bystander_exposure = pmin(1,prob_bystander_exposure)
  )

# df_bacteria <- filter_coexistence_condition(df_bacteria)



#### Run analysis equilibrium ####
cat("Computing coexistence equilibrium...")

# Compute equilibrium with analytical expressions 
eq_results <- compute_equilibrium(df_bacteria, population_size)

# # Check stability of equilibrium
# eq_results <- eq_results %>%
#   mutate(stability = pmap(., check_eq_stability)) %>%
#   unnest_wider(stability) 
#   
# 
# create_folder(file.path(getwd(),"files",folder_name))
save(eq_results, file = here::here("files",folder_name,"equilibrium_results_2.RData"))

cat("done and equilibrium results saved\n")


#### Simulation de un an sans vaccin et un an avec vaccin ####
# Simulate 1 year without vaccine for each parameter set at equilibrium
cat("Computing 1 year without vaccination...")

results_1y_wov <- eq_results %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov") %>%
  unnest_wider(res_1y_wov, names_sep = "_")

cat("done\n")

cat("Computing 1 year with vaccination...")
# Choose vaccination scenario
vaccine_scenarios_complete_df <- bind_rows(
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=0,"vfdr"=0,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0,name = "vftc"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vfds" = c(0.3,0.6,0.9)) %>% mutate("vfdr" = vfds, "vftcs"=0, "vftcr"=0,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0,name = "vfd"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vfis" = c(0.3,0.6,0.9)) %>% mutate("vfir" = vfis,"vftcs"=0, "vftcr"=0, "vfds"=0, "vfdr"=0, "vfrs"=0, "vfrr"=0, name = "vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfd"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfds"=0,"vfdr"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vfds" = c(0.3,0.6,0.9)) %>% mutate("vfdr" = vfds, "vfis"=vfds,"vfir"=vfds,"vftcs"=0, "vftcr"=0,"vfrs"=0, "vfrr"=0, name = "vfd_vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfrs"=0, "vfrr"=0, name = "vftc_vfd_vfi")
) %>% mutate(vaccine_id = row_number())

if(bacteria == "E_coli"){
  vaccine_scenarios_complete_df <- vaccine_scenarios_complete_df %>%
    bind_rows(tibble("Vperc" = c(0.7, 0.7), 
                     "vftcs" = c(0.7,0),"vfds" = c(0.7, 0), "vfis" = c(0.7, 0.7), "vfrs" = c(0,0),
                     "vftcr" = c(0.7,0),"vfdr" = c(0.7, 0), "vfir" = c(0.7, 0.7), "vfrr" = c(0,0),
                     "name" = c("vftc_vfd_vfi", "vfi"), 
                     "vaccine_id" =c(max(vaccine_scenarios_complete_df$vaccine_id)+1,max(vaccine_scenarios_complete_df$vaccine_id)+2) ))
}


if(transmission_by_infected){
  vaccine_scenarios_complete_df <- vaccine_scenarios_complete_df %>%
    mutate(vftis = vftcs, vftir = vftcr)
} else {
  vaccine_scenarios_complete_df <- vaccine_scenarios_complete_df %>%
    mutate(vftis = 0, vftir = 0)
}



# Add vaccine parameters to each parameter set
df_for_vaccine_simulations <- add_vaccine_parameters(results_1y_wov, vaccine_scenarios_complete_df)

# df_for_vaccine_simulations_test <- df_for_vaccine_simulations %>%
#   mutate(coexistence_condition_vaccine = f <= (1/(dpr*(1-vfdr)) + ar*(1-vfir) - etar * ar*(1-vfir) / (etar + gammaA_r)) /(1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas * as*(1-vfis) / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as*(1-vfis)/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar*(1-vfir)/(etar + gammaA_r))) %>%
#   mutate(non_disparition_condition_vaccine = (1-vftcs)*(betaC + betaI*as*(1-vfis)/(etas + lambdaA_I + gammaA_s + psiA + phiA_I)) >= 1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas * as*(1-vfis) / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)%>%
#   mutate(
#     eq_Sv = N/(betaC*(1-vftcs) + betaI*(1-vftcs)*as*(1-vfis)/His)*(1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas*as*(1-vfis)/His + phiA),
#     E1v = f*eq_Sv/N*(betaC*(1-vftcr) + betaI*(1-vftcr)*ar*(1-vfir)/(etar+gammaA_r))-1/(dpr*(1-vfdr))- ar*(1-vfir) + etar*ar*(1-vfir)/(etar+gammaA_r),
#     E2v = (betaI*(1-vftcr)*f*eq_Sv/N +etar)*(psiA + phiA_I)/(etar+gammaA_r)*as*(1-vfis)/His +phiA,
#     Ev = E1v/E2v,
#     Gv = 1+ ar*(1-vfir)/(etar+gammaA_r)-Ev*(1+as*(1-vfis)/His*(1+(psiA + phiA_I)/(etar+gammaA_r))),
#     eq_Crv = (N-eq_Sv)/Gv,
#     eq_Csv = -Ev*eq_Crv
#   )%>%
#   mutate(test1 = 1/(1-E),
#          test2 = 1/(1-Ev)) %>%
#   mutate(test = test2 - test1)

# Simulate one year with vaccine for each parameter set
results <- df_for_vaccine_simulations %>%
  apply_function_on_df(simulate_1y_with_vaccine, "res_1y_wv") %>%
  unnest_wider(res_1y_wv, names_sep = "_")


#### Results manipulation ####
# Compute resistance change
results_with_outputs <- compute_outputs(results, population_size)

results_with_outputs <- clean_vaccine_related_parameters(results_with_outputs)

# save results_with_outputs
save(results_with_outputs, file = here::here("files",folder_name,"results_with_outputs_2.RData"))

# load results_with_outputs
#load(here::here("files",folder_name,"results_with_outputs.RData"))

# Compute statistics
results_to_plot <- results_with_outputs %>%
  compute_statistics(starts_with(c("res_1y_wv_prevC","res_1y_wv_inccumI","prc_red_", "averted_")))%>%
  mutate(
    population = case_when(
      str_detect(metric_name, "non_vaccinated") ~ "Non-vaccinated",
      str_detect(metric_name, "_vaccinated") ~ "Vaccinated",
      TRUE ~ "Total population"
    )
  ) %>%
  mutate(population = factor(population, levels = c("Total population", "Vaccinated", "Non-vaccinated")))

# save results_to_plot
save(results_to_plot, file = here::here("files",folder_name,"results_to_plot_2.RData"))
cat("...done and results saved")
