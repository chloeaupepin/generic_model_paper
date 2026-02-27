#### Reproduction number ####


SCISsrV.reprod_nbrs <- function(param){
  with(as.list(c(param)), {
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
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
    Hcsv = 1/dps*1/(1-vfds)+lambdaA+phiA+as*(1-vfis)
    Hisnv = etas+lambdaA_I+gammaA_s+psiA+phiA_I
    Hisv = etas*1/(1-vfrs)+lambdaA_I+gammaA_s+psiA+phiA_I
    Hcrnv = 1/dpr+ar
    Hcrv = 1/dpr*1/(1-vfdr)+ar*(1-vfir)
    Hirnv = etar+gammaA_r
    Hirv = etar*1/(1-vfrr)+gammaA_r
    
    Pnv = 1-Vperc
    Pv = Vperc
    
    detsnv = Hcsnv*Hisnv-as*(etas+lambdaA_I*(1-eps))
    detsv = Hcsv*Hisv-as*(1-vfis)*(etas*1/(1-vfrs)+lambdaA_I*(1-eps))
    detrnv = Hcrnv*Hirnv-ar*etar
    detrv = Hcrv*Hirv-ar*(1-vfir)*etar*1/(1-vfrr)
    
    X_inv = matrix(c(Hisnv/detsnv, 0 , (etas+lambdaA_I*(1-eps))/detsnv, 0,
                   0, Hisv/detsv, 0, (etas*1/(1-vfrs)+lambdaA_I*(1-eps))/detsv,
                   as/detsnv, 0, Hcsnv/detsnv, 0,
                   0, as*(1-vfis)/detsv,0,Hcsv/detsv), nrow=4, byrow = TRUE)
    
    Y_inv = matrix(c(Hirnv/detrnv, 0 , etar/detrnv, 0,
                   0, Hirv/detrv, 0, etar*1/(1-vfrr)/detrv,
                   ar/detrnv, 0, Hcrnv/detrnv, 0,
                   0, ar*(1-vfir)/detrv,0,Hcrv/detrv), nrow=4, byrow = TRUE)
    
    Z = diag(x=c(phiA,phiA,psiA+phiA_I,psiA+phiA_I),nrow = 4,ncol=4)
    
    Y_invZX_inv = Y_inv %*% Z %*% X_inv
    
    V_inv <- - rbind(cbind(X_inv,matrix(0,nrow=4,ncol=4)), cbind(Y_invZX_inv, Y_inv))
    Fm <- matrix(c(betaC*Pnv, betaC*Pnv, betaI*Pnv, betaI*Pnv, 0, 0, 0, 0,
                   betaC*(1-vftcs)*Pv,betaC*(1-vftcs)*Pv,betaI*(1-vftis)*Pv, betaI*(1-vftis)*Pv, 0,0,0,0,
                   0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,
                   0,0,0,0,betaC*f*Pnv, betaC*f*Pnv, betaI*f*Pnv, betaI*f*Pnv,
                   0,0,0,0,betaC*f*(1-vftcr)*Pv,betaC*f*(1-vftcr)*Pv,betaI*f*(1-vftir)*Pv, betaI*f*(1-vftir)*Pv,
                   0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0), nrow = 8, byrow = TRUE)
    
    G = - Fm %*% V_inv
    
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


vaccine_coverage_threshold_for_R0 <- function(param){
  with(as.list(c(param)), {
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
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
    Hcsv = 1/dps*1/(1-vfds)+lambdaA+phiA+as*(1-vfis)
    Hisnv = etas+lambdaA_I+gammaA_s+psiA+phiA_I
    Hisv = etas*1/(1-vfrs)+lambdaA_I+gammaA_s+psiA+phiA_I
    Hcrnv = 1/dpr+ar
    Hcrv = 1/dpr*1/(1-vfdr)+ar*(1-vfir)
    Hirnv = etar+gammaA_r
    Hirv = etar*1/(1-vfrr)+gammaA_r
    
    detsnv = Hcsnv*Hisnv-as*(etas+lambdaA_I*(1-eps))
    detsv = Hcsv*Hisv-as*(1-vfis)*(etas*1/(1-vfrs)+lambdaA_I*(1-eps))
    detrnv = Hcrnv*Hirnv-ar*etar
    detrv = Hcrv*Hirv-ar*(1-vfir)*etar*1/(1-vfrr)
    
    k = (betaC*Hisnv + betaI*as)/detsnv
    l = (1-vftcs)*(betaC*Hisv + betaI*as*(1-vfis))/detsv
    # m = (1-vftcs)*(betaC*Hisnv + betaI*as)/detsnv * (betaC*Hisv + betaI*as*(1-vfis))/detsv
    # print(k)
    # print(l)
    # print(m)
    # 
    # print((k*l+k-l-m)^2 - 4*(1-k)*(m-k*l))
    # print((m-k*l))
    # 
    # Vperc1 = (-sqrt((k*l+k-l-m)^2 - 4*(1-k)*(m-k*l))-k*l-k+l+m)/(2*(m-k*l))
    # Vperc2 = (sqrt((k*l+k-l-m)^2 - 4*(1-k)*(m-k*l))-k*l-k+l+m)/(2*(m-k*l))
    
    Vperc = (1-k)/(l-k)
      
    # return(tibble("Vperc1"=Vperc1, "Vperc2"=Vperc2))
    return(Vperc)
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
  data <- compute_antibiotic_associated_flux(data, is.null(data$gammaA_s))
  
  
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
apply_function_on_df <- function(df, fun, col_name, map_fun = future_pmap) {
  df %>%
    mutate(!!col_name := map_fun(., function(...) {
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
    left_join(df_to_add, by = "vaccine_id")
  return(df)
}

# Compute outputs 
prc_red <- function(ref, new) {
  -(ref - new) / ref * 100
}

compute_outputs <- function(sim, population_size = 100000, factor = 100000){
  # Add cumulative incidence of selection to cumulative incidence of resistant infections
  results <- sim %>%
    mutate(res_1y_wov_inccumIrnv_sel = res_1y_wov_inccumIrnv + res_1y_wov_inccumSelectionOfResistantBySpecificAntibioForNonVaccinated + res_1y_wov_inccumSelectionOfResistantByBystanderAntibioForNonVaccinated,
           res_1y_wv_inccumIrnv_sel = res_1y_wv_inccumIrnv + res_1y_wv_inccumSelectionOfResistantBySpecificAntibioForNonVaccinated + res_1y_wv_inccumSelectionOfResistantByBystanderAntibioForNonVaccinated,
           res_1y_wv_inccumIrv_sel = res_1y_wv_inccumIrv + res_1y_wv_inccumSelectionOfResistantBySpecificAntibioForVaccinated + res_1y_wv_inccumSelectionOfResistantByBystanderAntibioForVaccinated)
    
  
  # Add cumulative incidence of total infections
  results <- results %>%
    mutate(res_1y_wov_inccumI = res_1y_wov_inccumIrnv_sel+res_1y_wov_inccumIsnv,
           res_1y_wv_inccumI_non_vaccinated = res_1y_wv_inccumIrnv_sel + res_1y_wv_inccumIsnv,
           res_1y_wv_inccumI_vaccinated = res_1y_wv_inccumIrv_sel + res_1y_wv_inccumIsv,
           res_1y_wv_inccumI = res_1y_wv_inccumI_non_vaccinated + res_1y_wv_inccumI_vaccinated)
  
  # Compute reduction change in cumulative incidences
  results<- results %>%
    mutate(prc_red_inccumIs = prc_red(res_1y_wov_inccumIsnv,res_1y_wv_inccumIsnv + res_1y_wv_inccumIsv),
           prc_red_inccumIs_non_vaccinated = prc_red(res_1y_wov_inccumIsnv*(1-Vperc),res_1y_wv_inccumIsnv ),
           prc_red_inccumIs_vaccinated = prc_red(res_1y_wov_inccumIsnv*Vperc, res_1y_wv_inccumIsv),
           prc_red_inccumIr = prc_red(res_1y_wov_inccumIrnv_sel,res_1y_wv_inccumIrnv_sel + res_1y_wv_inccumIrv_sel),
           prc_red_inccumIr_non_vaccinated = prc_red(res_1y_wov_inccumIrnv_sel*(1-Vperc),res_1y_wv_inccumIrnv_sel ),
           prc_red_inccumIr_vaccinated = prc_red(res_1y_wov_inccumIrnv_sel*Vperc, res_1y_wv_inccumIrv_sel),
           prc_red_inccumI = prc_red(res_1y_wov_inccumI, res_1y_wv_inccumI),
           prc_red_inccumI_non_vaccinated = prc_red(res_1y_wov_inccumI*(1-Vperc), res_1y_wv_inccumI_non_vaccinated),
           prc_red_inccumI_vaccinated = prc_red(res_1y_wov_inccumI*Vperc, res_1y_wv_inccumI_vaccinated)
           )
  
  # Compute proportion of resistant infections among infections in cumulative incidences
  results <- results %>%
    mutate(res_1y_wov_prop_inccumIr = (res_1y_wov_inccumIrnv_sel)/res_1y_wov_inccumI,
           res_1y_wv_prop_inccumIr = (res_1y_wv_inccumIrnv_sel + res_1y_wv_inccumIrv_sel)/res_1y_wv_inccumI,
           res_1y_wv_prop_inccumIr_non_vaccinated = (res_1y_wv_inccumIrnv_sel)/res_1y_wv_inccumI_non_vaccinated,
           res_1y_wv_prop_inccumIr_vaccinated = (res_1y_wv_inccumIrv_sel)/res_1y_wv_inccumI_vaccinated
    )
  
  
  # Compute proportion of resistant colonisation among colonisation in prevalences at the end of the year
  results <- results %>%
    mutate(res_1y_wv_prevC =( res_1y_wv_Crnv + res_1y_wv_Csnv + res_1y_wv_Crv  + res_1y_wv_Csv)/population_size) %>%
    mutate(res_1y_wov_prop_prevCIr = (res_1y_wov_Crnv+res_1y_wov_Irnv)/(res_1y_wov_Crnv + res_1y_wov_Irnv + res_1y_wov_Csnv + res_1y_wov_Isnv),
           res_1y_wv_prop_prevCIr = (res_1y_wv_Crnv+res_1y_wv_Irnv + res_1y_wv_Crv+res_1y_wv_Irv)/(res_1y_wv_Crnv + res_1y_wv_Irnv + res_1y_wv_Csnv + res_1y_wv_Isnv + res_1y_wv_Crv + res_1y_wv_Irv + res_1y_wv_Csv + res_1y_wv_Isv),
           res_1y_wv_prop_prevCIr_non_vaccinated = (res_1y_wv_Crnv+res_1y_wv_Irnv)/(res_1y_wv_Crnv + res_1y_wv_Irnv + res_1y_wv_Csnv + res_1y_wv_Isnv),
           res_1y_wv_prop_prevCIr_vaccinated = (res_1y_wv_Crv+res_1y_wv_Irv)/(res_1y_wv_Crv + res_1y_wv_Irv + res_1y_wv_Csv + res_1y_wv_Isv),
           
           res_1y_wov_prop_prevCr = (res_1y_wov_Crnv)/(res_1y_wov_Crnv + res_1y_wov_Csnv),
           res_1y_wv_prop_prevCr = (res_1y_wv_Crnv + res_1y_wv_Crv)/(res_1y_wv_Crnv + res_1y_wv_Csnv + res_1y_wv_Crv  + res_1y_wv_Csv ),
           res_1y_wv_prop_prevCr_non_vaccinated = (res_1y_wv_Crnv)/(res_1y_wv_Crnv + res_1y_wv_Csnv ),
           res_1y_wv_prop_prevCr_vaccinated = (res_1y_wv_Crv)/(res_1y_wv_Crv + res_1y_wv_Csv )
           )
  
  # Compute reduction change in proportions
  results<- results %>%
    mutate(prc_red_prop_prevCIr = prc_red(res_1y_wov_prop_prevCIr,res_1y_wv_prop_prevCIr ),
           prc_red_prop_prevCIr_non_vaccinated = prc_red(res_1y_wov_prop_prevCIr,res_1y_wv_prop_prevCIr_non_vaccinated ),
           prc_red_prop_prevCIr_vaccinated = prc_red(res_1y_wov_prop_prevCIr,res_1y_wv_prop_prevCIr_vaccinated ))%>%
    mutate(prc_red_prop_prevCr = prc_red(res_1y_wov_prop_prevCr,res_1y_wv_prop_prevCr ),
           prc_red_prop_prevCr_non_vaccinated = prc_red(res_1y_wov_prop_prevCr,res_1y_wv_prop_prevCr_non_vaccinated ),
           prc_red_prop_prevCr_vaccinated = prc_red(res_1y_wov_prop_prevCr,res_1y_wv_prop_prevCr_vaccinated ))%>%
    mutate(prc_red_prop_inccumIr = prc_red(res_1y_wov_prop_inccumIr,res_1y_wv_prop_inccumIr ),
           prc_red_prop_inccumIr_non_vaccinated = prc_red(res_1y_wov_prop_inccumIr,res_1y_wv_prop_inccumIr_non_vaccinated ),
           prc_red_prop_inccumIr_vaccinated = prc_red(res_1y_wov_prop_inccumIr,res_1y_wv_prop_inccumIr_vaccinated ))
  
  # Compute reduction change in prevalences
  results <- results %>%
    mutate(prc_red_prevCr = prc_red(res_1y_wov_Crnv, res_1y_wv_Crnv + res_1y_wv_Crv),
           prc_red_prevCr_non_vaccinated = prc_red(res_1y_wov_Crnv*(1-Vperc), res_1y_wv_Crnv),
           prc_red_prevCr_vaccinated = prc_red(res_1y_wov_Crnv*Vperc, res_1y_wv_Crv))%>%
    mutate(prc_red_prevC = prc_red(res_1y_wov_Crnv+res_1y_wov_Csnv,res_1y_wv_Crnv + res_1y_wv_Crv + res_1y_wv_Csnv + res_1y_wv_Csv),
           prc_red_prevC_non_vaccinated = prc_red((res_1y_wov_Crnv+res_1y_wov_Csnv)*(1-Vperc),res_1y_wv_Crnv +res_1y_wv_Csnv ),
           prc_red_prevC_vaccinated = prc_red((res_1y_wov_Crnv+res_1y_wov_Csnv)*Vperc, res_1y_wv_Crv + res_1y_wv_Csv))
  
  
  # Compute averted infections 
  results <- results %>%
    mutate(averted_inccumIs = res_1y_wov_inccumIsnv - res_1y_wv_inccumIsnv - res_1y_wv_inccumIsv,
           averted_inccumIr = res_1y_wov_inccumIrnv_sel - res_1y_wv_inccumIrnv_sel - res_1y_wv_inccumIrv_sel,
           averted_inccumI = averted_inccumIs + averted_inccumIr,
           averted_inccumIs_non_vaccinated = res_1y_wov_inccumIsnv*(1-Vperc) - res_1y_wv_inccumIsnv ,
           averted_inccumIr_non_vaccinated = res_1y_wov_inccumIrnv_sel*(1-Vperc) - res_1y_wv_inccumIrnv_sel,
           averted_inccumI_non_vaccinated = averted_inccumIs_non_vaccinated + averted_inccumIr_non_vaccinated,
           averted_inccumIs_vaccinated = res_1y_wov_inccumIsnv*Vperc - res_1y_wv_inccumIsv,
           averted_inccumIr_vaccinated = res_1y_wov_inccumIrnv_sel*Vperc - res_1y_wv_inccumIrv_sel,
           averted_inccumI_vaccinated = averted_inccumIs_vaccinated + averted_inccumIr_vaccinated)%>%
    mutate(across(starts_with("averted"), ~ .x * factor/population_size))
  
  return(results)
  
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

create_folder <- function(path){
  if(!dir.exists(path)){
    dir.create(path)
  }
}
