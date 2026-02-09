#### Reproduction number ####


SCISsrV.reprod_nbrs <- function(param){
  with(as.list(c(param)), {
    
    #vftcr = vftcs
    #vfir = vfis
    #vfdr = vfds
    #vfrr = vfrs
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    lambdaA_I = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
    if(!("gammaA_s" %in% names(param))){
      gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*(1-prob_minority_strain_when_infected)
    }
    if(!("psiA" %in% names(param))){
      psiA = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*prob_minority_strain_when_infected
    }
    if(!("gammaA_r" %in% names(param))){
      gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    }
    
    etas = (1-prob_specific_exposure) * 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    Hcsnv = 1/dps+lambdaA+phiA+as
    Hcsv = 1/dps*(1+vfds)+lambdaA+phiA+as*(1-vfis)
    Hisnv = etas+lambdaA_I+gammaA_s+psiA+phiA_I
    Hisv = etas*(1+vfrs)+lambdaA_I+gammaA_s+psiA+phiA_I
    Hcrnv = 1/dpr+ar
    Hcrv = 1/dpr*(1+vfdr)+ar*(1-vfir)
    Hirnv = etar+gammaA_r
    Hirv = etar*(1+vfrr)+gammaA_r
    
    #Nnv = Snv0+Csnv0+Crnv0+Isnv0 + Irnv0
    #Nv = Sv0+Csv0+Crv0 + Isv0 + Irv0
    Pnv = Nnv/(Nnv+Nv)
    Pv = Nv/(Nnv+Nv)
    
    detsnv = Hcsnv*Hisnv-as*(etas+lambdaA_I*(1-eps))
    detsv = Hcsv*Hisv-as*(1-vfis)*(etas+lambdaA_I*(1-eps))
    detrnv = Hcrnv*Hirnv-ar*etar
    detrv = Hcrv*Hirv-ar*(1-vfir)*etar
    
    X_1 = matrix(c(Hisnv/detsnv, 0 , (etas+lambdaA_I*(1-eps))/detsnv, 0,
                   0, Hisv/detsv, 0, (etas+lambdaA_I*(1-eps))/detsv,
                   as/detsnv, 0, Hcsnv/detsnv, 0,
                   0, as*(1-vfis)/detsv,0,Hcsv/detsv), nrow=4, byrow = TRUE)
    
    Y_1 = matrix(c(Hirnv/detrnv, 0 , etar/detrnv, 0,
                   0, Hirv/detrv, 0, etar/detrv,
                   ar/detrnv, 0, Hcrnv/detrnv, 0,
                   0, ar*(1-vfir)/detrv,0,Hcrv/detrv), nrow=4, byrow = TRUE)
    
    Z = diag(x=c(phiA,phiA,psiA+phiA_I,psiA+phiA_I),nrow = 4,ncol=4)
    
    Y_1ZX_1 = Y_1 %*% Z %*% X_1
    
    Vm <- rbind(cbind(X_1,matrix(0,nrow=4,ncol=4)), cbind(Y_1ZX_1, Y_1))
    Fm <- matrix(c(betaC*Pnv, betaC*Pnv, betaI*Pnv, betaI*Pnv, 0, 0, 0, 0,
                   betaC*(1-vftcs)*Pv,betaC*(1-vftcs)*Pv,betaI*(1-vftis)*Pv, betaI*(1-vftis)*Pv, 0,0,0,0,
                   0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,
                   0,0,0,0,betaC*f*Pnv, betaC*f*Pnv, betaI*f*Pnv, betaI*f*Pnv,
                   0,0,0,0,betaC*f*(1-vftcr)*Pv,betaC*f*(1-vftcr)*Pv,betaI*f*(1-vftir)*Pv, betaI*f*(1-vftir)*Pv,
                   0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0), nrow = 8, byrow = TRUE)
    
    G = Fm %*% Vm
    
    rownames(G) <- c("de type Csnv","de type Csv","de type Isnv","de type Isv",
                     "de type Crnv","de type Crv","de type Irnv","de type Irv")
    
    colnames(G) <- c("par un individu de type Csnv","par un individu de type Csv","par un individu de type Isnv","par un individu de type Isv",
                     "par un individu de type Crnv","par un individu de type Crv","par un individu de type Irnv","par un individu de type Irv")
    
    df <- as.data.frame(G, check.names = FALSE)
    
    df <- cbind("Infections secondaires" = rownames(df),df)
    rownames(df) <- NULL
    
    vaps = 1/2*(G[1,1]+ G[2,2]+sqrt(G[1,1]^2+G[2,2]^2-2*G[1,1]*G[2,2]+4*G[2,1]*G[1,2]))
    vapr = 1/2*(G[5,5]+ G[6,6]+sqrt(G[5,5]^2+G[6,6]^2-2*G[5,5]*G[6,6]+4*G[6,5]*G[5,6]))
    
    R0 <- max(vaps,vapr)
    return(tibble("R0s"= vaps, "R0r" = vapr))
    #return(list(df,R0))
  })
  
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
    left_join(vaccine_scenarios_complete_df, by = "vaccine_id")
  return(df)
}

