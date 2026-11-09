library(purrr)
library(tidyr)
library(dplyr)
library(deSolve)
library(furrr)
library(stringr)

source(here::here("scripts", "model.R"))

n_cores = parallel::detectCores() - 3
plan(multisession, workers = n_cores)

set.seed(123)  # global seed for reproducibility



#### Helpers ####

# Apply a function on dataframe rows
apply_function_on_df <- function(df, fun, col_name, map_fun = future_pmap) {
  df %>%
    mutate(!!col_name := map_fun(., function(...) {
      params <- list(...)
      fun(params)
    }))
}


create_folder <- function(path){
  if(!dir.exists(path)){
    dir.create(path)
  }
}


#### Compute antibiotic associated flux ####

compute_antibiotic_associated_flux <- function(data, with_specific_flux=T){
  if(with_specific_flux){
    data <- data %>%
      mutate(gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*(1-prob_minority_strain_when_infected),
             psiA = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*prob_minority_strain_when_infected,
             gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r,
             )
  }
  data %>%
    mutate(lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised),
           phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised,
           
           lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected),
           phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected,
           
           etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s,
           etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r)
}

#### Fix resistant parameters equal to sensitive ####

set_resistant_and_sensitive_param_equal <- function(data, res_param_names){
  data_with_fixed_res_param <- data
  if("dpr" %in% res_param_names ){
    data_with_fixed_res_param <- data_with_fixed_res_param %>% 
      mutate(dpr = dps)
  }
  if("ar" %in% res_param_names ){
    data_with_fixed_res_param <- data_with_fixed_res_param %>% 
      mutate(ar = as)
  }
  if("time_until_recovery_without_ATB_r" %in% res_param_names ){
    data_with_fixed_res_param <- data_with_fixed_res_param %>% 
      mutate(time_until_recovery_without_ATB_r = time_until_recovery_without_ATB_s)
  }
  if("vftcr" %in% res_param_names ){
    data_with_fixed_res_param <- data_with_fixed_res_param %>% 
      mutate(vftcr = vftcs)
  }
  if("vftir" %in% res_param_names ){
    data_with_fixed_res_param <- data_with_fixed_res_param %>% 
      mutate(vftir = vftis)
  }
  if("vfdr" %in% res_param_names ){
    data_with_fixed_res_param <- data_with_fixed_res_param %>% 
      mutate(vfdr = vfds)
  }
  if("vfir" %in% res_param_names ){
    data_with_fixed_res_param <- data_with_fixed_res_param %>% 
      mutate(vfir = vfis)
  }
  if("vfrr" %in% res_param_names ){
    data_with_fixed_res_param <- data_with_fixed_res_param %>% 
      mutate(vfrr = vfrs)
  }
  
  return(data_with_fixed_res_param)
}



#### Analytical conditions for coexistence ####

compute_analytical_conditions <- function(data){
  data <- compute_antibiotic_associated_flux(data, is.null(data$gammaA_s))
  
  
  data %>%
    mutate(coexistence_condition = f <= (1/dpr + ar - etar * ar / (etar + gammaA_r)) /(1/dps + lambdaA + as - etas * as / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar/(etar + gammaA_r))) %>%
    mutate(non_disparition_condition = betaC + betaI*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I) >= 1/dps + lambdaA + as - etas * as / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)
}

filter_coexistence_condition <- function(data){
  compute_analytical_conditions(data) %>%
    filter(coexistence_condition & non_disparition_condition)
}


#### Compute equilibrium ####
compute_equilibrium <- function(data,population_size = 100000){
  data <- compute_antibiotic_associated_flux(data, is.null(data$gammaA_s))
  
  
  data %>%
    mutate(
      N = population_size,
      His = etas + lambdaA_I +  gammaA_s+ psiA + phiA_I,
      eq_Snv = N/(betaC + betaI*as/His)*(1/dps + lambdaA + as - etas*as/His + phiA),
      E1 = f*eq_Snv/N*(betaC + betaI*ar/(etar+gammaA_r))-1/dpr - ar + etar*ar/(etar+gammaA_r),
      E2 = (betaI*f*eq_Snv/N +etar)*(psiA + phiA_I)/(etar+gammaA_r)*as/His +phiA,
      E = E1/E2,
      G = 1+ ar/(etar+gammaA_r)-E*(1+as/His*(1+(psiA + phiA_I)/(etar+gammaA_r))),
      eq_Crnv = (N-eq_Snv)/G,
      eq_Csnv = -E*eq_Crnv,
      eq_Isnv = as/His*eq_Csnv,
      eq_Irnv = ar/(etar+gammaA_r)*eq_Crnv + (psiA + phiA_I)/(etar+gammaA_r)*eq_Isnv
    )
}

check_eq_stability <- function(...){
  
  param <- list(...)
  
  with(param, {
    
    J = matrix(c(-(betaC*eq_Csnv + betaI*eq_Isnv)/N - f*(betaC*eq_Crnv + betaI*eq_Irnv)/N,
                 -betaC*eq_Snv/N + 1/dps +lambdaA,
                 -f*betaC*eq_Snv/N + 1/dpr,
                 -betaI*eq_Snv/N + lambdaA_I + gammaA_s,
                 -f*betaI*eq_Snv/N + gammaA_r,
                 
                 (betaC*eq_Csnv + betaI*eq_Isnv)/N,
                 betaC*eq_Snv/N - 1/dps - lambdaA - as - phiA,
                 0,
                 betaI*eq_Snv/N + etas,
                 0,
                 
                 f*(betaC*eq_Crnv + betaI*eq_Irnv)/N,
                 phiA,
                 f*betaC*eq_Snv/N - 1/dpr - ar,
                 0,
                 f*betaI*eq_Snv/N + etar,
                 
                 0,
                 as,
                 0,
                 -etas - lambdaA_I - gammaA_s - psiA - phiA_I,
                 0,
                 
                 0,
                 0,
                 ar,
                 psiA + phiA_I,
                 -etar - gammaA_r
                 ), nrow=5, byrow = TRUE)
    
    
    J = matrix(c(betaC*eq_Snv/N - (betaC*eq_Csnv + betaI*eq_Isnv)/N - 1/dps - lambdaA - as - phiA,
                 - (betaC*eq_Csnv + betaI*eq_Isnv)/N,
                 betaI*eq_Snv/N + etas- (betaC*eq_Csnv + betaI*eq_Isnv)/N,
                 - (betaC*eq_Csnv + betaI*eq_Isnv)/N,
                 
                 
                 phiA - f*(betaC*eq_Crnv + betaI*eq_Irnv)/N,
                 f*betaC*eq_Snv/N - f*(betaC*eq_Crnv + betaI*eq_Irnv)/N - 1/dpr - ar,
                 - f*(betaC*eq_Crnv + betaI*eq_Irnv)/N,
                 f*betaI*eq_Snv/N - f*(betaC*eq_Crnv + betaI*eq_Irnv)/N + etar,
                 
                 
                 as,
                 0,
                 -etas - lambdaA_I - gammaA_s - psiA - phiA_I,
                 0,
                 
                 
                 0,
                 ar,
                 psiA + phiA_I,
                 -etar - gammaA_r
    ), nrow=4, byrow = TRUE)
    
    eig = eigen(J)$values
    eig_neg = all(Re(eig) < 0)
    # print(eig)
    # print(eig_neg)
    
    #return(tibble("eig_neg"=eig_neg))
    #return(tibble("J" = as_tibble(J), "eig"= as_tibble(eig)) )
    return(tibble("eig_max"=max(eig)))
  })
}


#### Simulations ####

solve_model_equations <- function(y0, times, params){
  ode(
    y = y0,
    times = times,
    func = SCISsrV.model,
    parms = params
  ) %>% as.data.frame()
}

add_zero_init_cumulative_incidences <- function(init){
  c(init, c(inccumIsnv=0,inccumIsrnv=0, inccumIrnv=0, inccumIsv=0,inccumIsrv=0, inccumIrv=0, inccumSelectionOfResistantForNonVaccinated=0, inccumSelectionOfResistantForVaccinated=0))
}

arbitrary_initial_state <- function(){
  add_zero_init_cumulative_incidences(c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
                                        Sv = 0,Csv = 0,Csrv =0,Crv = 0,
                                        Isnv = 0,Isrnv = 0,Irnv = 0,
                                        Isv = 0,Isrv = 0,Irv = 0))
}

vaccine_params_zero <- function() {
  list(
    vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
    vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0
  )
}

# Function to run model until equilibrium without vaccine
simulate_eq_without_vaccine <- function(params) {
  out <- ode(
    y = arbitrary_initial_state(),
    times = seq(from=0,to=1000000,by=100),
    func = SCISsrV.model_3,
    parms = c(params, vaccine_params_zero())  
  ) %>% as.data.frame()
  
  out %>%
    slice_tail(n=1) %>% 
    select(Snv, Csnv,Csrnv, Crnv, Isnv,Isrnv, Irnv)
}

# Function to run model one year without vaccine

simulate_1y_without_vaccine <- function(params) {
  
  out <- ode(
    y = add_zero_init_cumulative_incidences(
      c(Snv = params$eq_Snv,Csnv = params$eq_Csnv,Csrnv = params$eq_Csrnv,Crnv = params$eq_Crnv,
          Sv = 0,Csv = 0, Csrv = 0, Crv =0, 
          Isnv = params$eq_Isnv ,Isrnv = params$eq_Isrnv ,Irnv = params$eq_Irnv , Isv = 0,Isrv = 0, Irv =0)),
    times = seq(from=0,to=370,by=10),
    func = SCISsrV.model_3,
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
    y = add_zero_init_cumulative_incidences(
      c(Snv = params$eq_Snv * (1-params$Vperc),Csnv = params$eq_Csnv * (1-params$Vperc),Csrnv = params$eq_Csrnv * (1-params$Vperc),Crnv = params$eq_Crnv * (1-params$Vperc),
          Sv = params$eq_Snv * (params$Vperc),Csv = params$eq_Csnv * (params$Vperc),Csrv = params$eq_Csrnv * (params$Vperc), Crv =params$eq_Crnv * (params$Vperc), 
          Isnv = params$eq_Isnv * (1-params$Vperc),Isrnv = params$eq_Isrnv * (1-params$Vperc),Irnv = params$eq_Irnv * (1-params$Vperc), Isv = params$eq_Isnv * (params$Vperc),Isrv = params$eq_Isrnv * (params$Vperc), Irv =params$eq_Irnv * (params$Vperc)
          )),
    times = seq(from=0,to=370,by=10),
    func = SCISsrV.model_3,
    parms = c(params)  
  ) %>% as.data.frame()
  
  out %>%
    slice_tail(n=1) %>% 
    select(Snv, Csnv,Csrnv, Crnv, Isnv,Isrnv, Irnv, Sv, Csv,Csrv, Crv, Isv,Isrv, Irv,
           inccumIsnv, inccumIsrnv, inccumIrnv, inccumIsv,inccumIsrv, inccumIrv,
           inccumSelectionOfResistantForNonVaccinated,inccumSelectionOfResistantForVaccinated)
}



#### Modify parameters ####

add_variability <- function(ref_values, params_to_change, n = 10, variation = 0.05) {
  # Create a matrix of random factors for the parameters to change
  random_factors <- matrix(runif(length(params_to_change) * n, 
                                 1 - variation, 1 + variation),
                           nrow = n, ncol = length(params_to_change))
  
  # Repeat the reference values n times and apply the random factors
  random_rows <- ref_values[rep(1, n), ]  # repeat the reference row n times
  random_rows[params_to_change] <- random_rows[params_to_change] * random_factors
  
  # Add bacteria IDs
  random_rows <- random_rows %>% mutate(bacteria_id = 1:n)
  
  # Add reference row with bacteria_id = 0 and put ID first
  df <- bind_rows(
    ref_values %>% mutate(bacteria_id = 0),
    random_rows
  ) %>% select(bacteria_id, everything()) %>%
    as_tibble()  # remove any weird row.names causing 1.1, 1.2 etc.
  
  return(df)
}

add_vaccine_parameters <- function(df_to_which_add,df_to_add){
  df <- df_to_which_add %>%
    slice(rep(1:n(), each = nrow(df_to_add))) %>%
    mutate(vaccine_id = rep(1:nrow(df_to_add), times = nrow(df_to_which_add))) %>%
    left_join(df_to_add, by = "vaccine_id")
  return(df)
}

clean_vaccine_related_parameters <- function(sim){
  # Add one column with each varying param
  results <- sim %>%
    mutate(varying_param = case_when(
      vftcs != 0 ~ vftcs,
      vfds != 0 ~ vfds,
      vfis != 0 ~ vfis,
      TRUE ~ 0
    )) %>% 
    filter(varying_param != 0) %>%
    mutate(name = factor(name, levels = c("vftc","vfd","vfi","vftc_vfd","vftc_vfi","vfd_vfi","vftc_vfd_vfi"))) %>%
    mutate(varying_param = factor(varying_param))
  
  # Rename vaccine names
  results <- results %>%
    mutate(name_renamed= case_when(
      name == "vftc" ~ "va",
      name == "vfd" ~ "vd",
      name == "vfi" ~ "vi",
      name == "vftc_vfd" ~ "va_vd",
      name == "vftc_vfi" ~ "va_vi",
      name == "vfd_vfi" ~ "vd_vi",
      name == "vftc_vfd_vfi" ~ "va_vd_vi"
    )) %>%
    mutate(name_renamed = factor(name_renamed, levels = c("va","vd","vi","va_vd","va_vi","vd_vi","va_vd_vi")))
  
  return(results)
}

#### Compute outputs ####

prc_red <- function(ref, new) {
  -(ref - new) / ref * 100
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

compute_statistics <- function(results, chosen_cols){
  results_with_stats <- results %>%
    pivot_longer(cols = all_of(chosen_cols), names_to = "metric_name", values_to = "metric_value")%>%
    group_by(name_renamed, Vperc, varying_param, metric_name) %>%
    summarise(median = median(metric_value),
              q025 = quantile(metric_value, 0.025),
              q975 = quantile(metric_value, 0.975),
              .groups = "drop")
  
  return(results_with_stats)
}


#### Define settings ####

## Population size
population_size = 100000

## Choose bacteria 
# Comment or uncomment 4 line codes below and rerun everything for each bacteria 
# S_aureus
bacteria = "S_aureus"
folder_name = "S_aureus"
file_name  = "S_aureus_params_eu_with_J01CR.csv"
transmission_by_infected = FALSE

# E_coli
bacteria = "E_coli"
folder_name = "E_coli"
file_name = "E_coli_params.csv"
transmission_by_infected = TRUE # if transmission is allowed, it equals transmission by colonized individuals

## Choose number of bacteria to generate
n <- 1


#### Define bacteria ####
Bacteria_params = read.csv(here::here("files",file_name)) 

# cols_to_noise include bacteria parameters and antibiotic parameters
# the relative fitness f isn't varied due to too much sensitivity of resistant proportion
if(bacteria == "S_aureus"){
  Bacteria_params$f = 0.94
  Bacteria_params$sigmaC1 = 1
  Bacteria_params$sigmaI1 = 0
  
  cols_to_noise <- c(	
    "betaC","as","time_until_recovery_without_ATB_s","dps",
    "prob_bystander_exposure",
    "time_until_decolo_by_bystander_ATB",
    "prob_minority_strain_when_colonised",
    "prob_specific_exposure",
    "time_until_decolo_by_specific_ATB",
    "prob_minority_strain_when_infected",
    "prob_specific_exposure_r",
    "time_until_decolo_by_specific_ATB_r","f"
  )
} else {
  
  Bacteria_params$sigmaC1 = 1.1e-1
  Bacteria_params$sigmaI1 = 1.1e-1
  
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

df_bacteria <- add_variability(Bacteria_params, cols_to_noise, n ,0.05) 

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


#Keep only first bacteria
df_bacteria <- df_bacteria %>% filter(bacteria_id == 0)

#### Run analysis equilibrium ####

# Compute equilibrium with analytical expressions 
eq_results <- df_bacteria %>%
  apply_function_on_df(simulate_eq_without_vaccine,"eq")%>%
  unnest_wider(eq, names_sep = "_")

eq_results_all <- eq_results

eq_results <- eq_results_all %>%
  mutate(prop_res = eq_Crnv/(eq_Csnv + eq_Csrnv + eq_Crnv)) 

if(bacteria =="S. aureus"){
  eq_results <- eq_results %>%
  filter(prop_res > 0.08)
}

save(eq_results, file = here::here("files",folder_name,"equilibrium_results_Csr_nequal1.RData"))


#### Simulation de un an sans vaccin et un an avec vaccin ####
# Simulate 1 year without vaccine for each parameter set at equilibrium

results_1y_wov <- eq_results %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov") %>%
  unnest_wider(res_1y_wov, names_sep = "_")

results_1y_wov_all <- results_1y_wov

results_1y_wov <- results_1y_wov_all %>%
  mutate(diff = eq_Snv - res_1y_wov_Snv + eq_Csnv - res_1y_wov_Csnv  + eq_Csrnv - res_1y_wov_Csrnv + eq_Crnv - res_1y_wov_Crnv + eq_Isnv - res_1y_wov_Isnv  + eq_Isrnv - res_1y_wov_Isrnv + eq_Irnv - res_1y_wov_Irnv) %>%
  filter(diff < 10-9)
  # filter(eq_Snv == res_1y_wov_Snv, 
  #        eq_Csnv == res_1y_wov_Csnv,eq_Csrnv == res_1y_wov_Csrnv, eq_Crnv == res_1y_wov_Crnv,
  #        eq_Isnv == res_1y_wov_Isnv,eq_Isrnv == res_1y_wov_Isrnv, eq_Irnv == res_1y_wov_Irnv)

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

# Simulate one year with vaccine for each parameter set
results <- df_for_vaccine_simulations %>%
  apply_function_on_df(simulate_1y_with_vaccine, "res_1y_wv") %>%
  unnest_wider(res_1y_wv, names_sep = "_")


#### Results manipulation ####
# Compute resistance change
results_with_outputs <- compute_outputs(results, population_size)

results_with_outputs <- clean_vaccine_related_parameters(results_with_outputs)


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

save(results_to_plot, file = here::here("files",folder_name,"results_to_plot_Csr_nequal1.RData"))

