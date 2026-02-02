#### Reproduction number ####


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
           
           etas = (1-prob_specific_exposure) * 1/time_until_recovery_without_ATB_s,
           etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r)
}


#### Analytical conditions for coexistence ####

compute_analytical_conditions <- function(data){
  if(is.null(data$gammaA_s)){
    data <- compute_antibiotic_associated_flux(data)
  }
  data %>%
    mutate(coexistence_condition = f <= (1/dpr + ar - etar * ar / (etar + gammaA_r)) /(1/dps + lambdaA + as - etas * as / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar/(etar + gammaA_r))) %>%
    mutate(non_disparition_condition = betaC + betaI*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I) >= 1/dps + lambdaA + as - etas * as / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)
}

filter_coexistence_condition <- function(data){
  compute_analytical_conditions(data) %>%
    filter(coexistence_condition & non_disparition_condition)
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

#### Compute equilibrium ####
compute_equilibrium <- function(data,population_size = 100000){
  if(is.null(data$gammaA_s)){
    data <- compute_antibiotic_associated_flux(data)
  }
  
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
    )%>%
    select(-N, -His, -E, -E1, -E2, -G) 
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
  c(init, c(inccumCsnv=0, inccumCrnv=0, inccumCsv=0, inccumCrv=0,
            inccumIsnv=0, inccumIrnv=0, inccumIsv=0, inccumIrv=0,
            
            inccumSensitiveNaturalRecoveryForNonVaccinated=0,inccumResistantNaturalRecoveryForNonVaccinated=0,
            inccumSensitiveNaturalRecoveryForVaccinated=0,inccumResistantNaturalRecoveryForVaccinated=0,
            
            inccumSensitiveSpecificAntibioRecoveryForNonVaccinated=0,inccumSensitiveSpecificAntibioRecoveryForVaccinated=0,
            inccumSensitiveBystanderAntibioRecoveryForNonVaccinated=0,inccumSensitiveBystanderAntibioRecoveryForVaccinated=0,
            
            inccumSelectionOfResistantBySpecificAntibioForNonVaccinated=0,inccumSelectionOfResistantBySpecificAntibioForVaccinated=0,
            inccumSelectionOfResistantByBystanderAntibioForNonVaccinated=0,inccumSelectionOfResistantByBystanderAntibioForVaccinated=0,
            
            inccumResistantSpecificAntibioRecoveryForNonVaccinated=0,inccumResistantSpecificAntibioRecoveryForVaccinated=0,
            
            inccumSuperInfectionByResistantForNonVaccinated=0,inccumSuperInfectionBySensitiveForNonVaccinated=0,
            inccumSuperInfectionByResistantForVaccinated=0,inccumSuperInfectionBySensitiveForVaccinated=0))
}

arbitrary_initial_state <- function(){
  add_zero_init_cumulative_incidences(c(Snv = 80000,Csnv = 10000,Crnv = 10000,Sv = 0,Csv = 0, Crv =0, Isnv = 0,Irnv = 0, Isv = 0, Irv =0))
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
    times = seq(from=0,to=80000,by=100),
    func = SCISsrV.model,
    parms = c(params, vaccine_params_zero())  
  ) %>% as.data.frame()
  
  out %>%
    slice_tail(n=1) %>% 
    select(Snv, Csnv, Crnv, Isnv, Irnv)
}

# Function to run model one year without vaccine

simulate_1y_without_vaccine <- function(params) {
  
  out <- ode(
    y = add_zero_init_cumulative_incidences(
      c(Snv = params$eq_Snv,Csnv = params$eq_Csnv,Crnv = params$eq_Crnv,
          Sv = 0,Csv = 0, Crv =0, 
          Isnv = params$eq_Isnv ,Irnv = params$eq_Irnv , Isv = 0, Irv =0)),
    times = seq(from=0,to=370,by=10),
    func = SCISsrV.model,
    parms = c(params, vaccine_params_zero())  
  ) %>% as.data.frame()
  
  out %>%
    slice_tail(n=1) %>% 
    select(Snv, Csnv, Crnv, Isnv, Irnv,
           inccumIsnv, inccumIrnv, inccumIsv, inccumIrv,
           inccumSelectionOfResistantBySpecificAntibioForNonVaccinated,
           inccumSelectionOfResistantByBystanderAntibioForNonVaccinated)
}

# Function to run model one year with vaccine

simulate_1y_with_vaccine <- function(params) {
  out <- ode(
    y = add_zero_init_cumulative_incidences(
      c(Snv = params$eq_Snv * (1-params$Vperc),Csnv = params$eq_Csnv * (1-params$Vperc),Crnv = params$eq_Crnv * (1-params$Vperc),
          Sv = params$eq_Snv * (params$Vperc),Csv = params$eq_Csnv * (params$Vperc), Crv =params$eq_Crnv * (params$Vperc), 
          Isnv = params$eq_Isnv * (1-params$Vperc),Irnv = params$eq_Irnv * (1-params$Vperc), Isv = params$eq_Isnv * (params$Vperc), Irv =params$eq_Irnv * (params$Vperc)
          )),
    times = seq(from=0,to=370,by=10),
    func = SCISsrV.model,
    parms = c(params)  
  ) %>% as.data.frame()
  
  out %>%
    slice_tail(n=1) %>% 
    select(Snv, Csnv, Crnv, Isnv, Irnv, Sv, Csv, Crv, Isv, Irv,
           inccumIsnv, inccumIrnv, inccumIsv, inccumIrv,
           inccumSelectionOfResistantBySpecificAntibioForNonVaccinated,
           inccumSelectionOfResistantByBystanderAntibioForNonVaccinated,
           inccumSelectionOfResistantBySpecificAntibioForVaccinated,
           inccumSelectionOfResistantByBystanderAntibioForVaccinated)
}

check_eq_without_vaccine <- function(params){
  with(as.list(c(params)), {
    Sv = 0
    Csv = 0
    Crv = 0
    Isv = 0
    Irv = 0
    Snv = eq_Snv
    Csnv = eq_Csnv
    Crnv = eq_Crnv
    Isnv = eq_Isnv
    Irnv = eq_Irnv
    
    N=Snv+Csnv+Crnv+Isnv+Irnv+Sv+Csv+Crv+Isv+Irv
    
    Bs = betaC*(Csnv+Csv)/N  + betaI*(Isnv + Isv)/N
    Br = betaC*f*(Crnv+Crv)/N + betaI*f*(Irnv+Irv)/N
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
    gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*(1-prob_minority_strain_when_infected)
    psiA = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*prob_minority_strain_when_infected
    
    gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    
    etas = (1-prob_specific_exposure) * 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    dSnv <- - Bs*Snv - Br*Snv + 1/dps*Csnv+1/dpr*Crnv +lambdaA*Csnv + gammaA_r*Irnv + lambdaA_I*eps*Isnv + gammaA_s*Isnv
    dCsnv <- Bs*Snv - 1/dps*Csnv - lambdaA*Csnv - phiA*Csnv - Br*Csnv*thetasr + Bs*Crnv*thetars -as*Csnv + etas*Isnv + lambdaA_I*(1-eps)*Isnv
    dIsnv <- as*Csnv - etas*Isnv - lambdaA_I*(1-eps)*Isnv - lambdaA_I*eps*Isnv - gammaA_s*Isnv - psiA*Isnv - phiA_I*Isnv
    dCrnv <- Br*Snv - 1/dpr*Crnv + phiA*Csnv + Br*Csnv*thetasr - Bs*Crnv*thetars - ar*Crnv + etar*Irnv
    dIrnv <- ar*Crnv - etar*Irnv + psiA*Isnv + phiA_I*Isnv -gammaA_r*Irnv 
    
    res<-sqrt(dSnv^2 + dCsnv^2 + dCrnv^2 + dIsnv^2+ dIrnv^2)
    
    res
  })
}

# Apply a function on dataframe rows
apply_function_on_df <- function(df, fun, col_name) {
  df %>%
    mutate(!!col_name := future_pmap(., function(...) {
      params <- list(...)
      fun(params)
    }))
}

