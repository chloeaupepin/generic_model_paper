
library(dplyr)
library(deSolve)
library(ggplot2)
library(tidyr)

SCISsrV.model_2 <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    # Total population size
    N=Snv+Csnv+Csrnv+Crnv+Isnv+Isrnv+Irnv+Sv+Csv+Csrv+Crv+Isv+Isrv+Irv
    
    ## Transmission coefficients
    Bs = betaC*(Csnv+Csv+Csrnv+Csrv)/N  + betaI*(Isnv+Isv+Isrnv+Isrv)/N
    Br = betaC*f*(Crnv+Crv)/N + betaI*f*(Irnv+Irv)/N
    Bsv = betaC*(1-vftcs)*(Csnv+Csv+Csrnv+Csrv)/N  + betaI*(1-vftis)*(Isnv+Isv+Isrnv+Isrv)/N
    Brv = betaC*(1-vftcr)*f*(Crnv+Crv)/N + betaI*(1-vftir)*f*(Irnv+Irv)/N
    
    
    # Bystander antibiotic exposure on colonised individuals
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    
    # Bystander antibiotic exposure on infected individuals
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    # Specific antibiotic exposure when infected
    if(!("gammaA_s" %in% names(param))){
      gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB
    }
    if(!("gammaA_r" %in% names(param))){
      gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    }
    
    # Recovery rates (when antibiotics aren't taken)
    etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    # Sigma 2 
    sigmaC2 = lambdaA
    sigmaI2 = lambdaA_I + gammaA_s
    
    # Equations of non-vaccinated compartments
    dSnv <- - Bs*Snv - Br*Snv + 1/dps*Csnv+1/dps*Csrnv+1/dpr*Crnv +lambdaA*Csnv + gammaA_r*Irnv + lambdaA_I*Isnv + gammaA_s*Isnv
    dCsnv <- Bs*Snv - 1/dps*Csnv - lambdaA*Csnv  -as*Csnv + etas*Isnv - sigmaC1*Csnv
    dCsrnv <- sigmaC1*Csnv - 1/dps*Csrnv - as*Csrnv - sigmaC2*Csrnv + etas*Isrnv
    dCrnv <- Br*Snv - 1/dpr*Crnv +  sigmaC2*Csrnv - ar*Crnv + etar*Irnv
    dIsnv <- as*Csnv - etas*Isnv - lambdaA_I*Isnv - gammaA_s*Isnv - sigmaI1*Isnv
    dIsrnv <- as*Csrnv - etas*Isrnv + sigmaI1*Isnv - sigmaI2*Isrnv
    dIrnv <- ar*Crnv - etar*Irnv -gammaA_r*Irnv + sigmaI2*Isrnv
    
    # Equations of vaccinated compartments
    dSv <- - Bsv*Sv - Brv*Sv + 1/dps*1/(1-vfds)*(Csv+Csrv)+1/dpr*1/(1-vfdr)*Crv + lambdaA*Csv + gammaA_r*Irv + lambdaA_I*Isv + gammaA_s*Isv
    dCsv <- Bsv*Sv - 1/dps*1/(1-vfds)*Csv - lambdaA*Csv - as*(1-vfis)*Csv + etas*1/(1-vfrs)*Isv - sigmaC1*Csv
    dCsrv <- sigmaC1*Csv - 1/dps*1/(1-vfds)*Csrv - as*(1-vfis)*Csrv - sigmaC2*Csrv + etas*1/(1-vfrs)*Isrv
    dCrv <- Brv*Sv  - 1/dpr*1/(1-vfdr)*Crv- ar*(1-vfir)*Crv + etar*1/(1-vfrr)*Irv + sigmaC2*Csrv
    dIsv <- as*(1-vfis)*Csv - etas*1/(1-vfrs)*Isv   - lambdaA_I*Isv - gammaA_s*Isv - sigmaI1*Isv
    dIsrv <- as*(1-vfis)*Csrv - etas*1/(1-vfrs)*Isrv + sigmaI1*Isv - sigmaI2*Isrv
    dIrv <- ar*(1-vfir)*Crv - etar*1/(1-vfrr)*Irv - gammaA_r*Irv + sigmaI2*Isrv

    
    res<-c(dSnv, dCsnv,dCsrnv, dCrnv,
           dSv, dCsv,dCsrv, dCrv, 
           dIsnv,dIsrnv, dIrnv, 
           dIsv,dIsrv, dIrv
           )
    
    list(res)
    
  })
  
}

get_equilibrium <- function(params) {
  
  with(as.list(c(params)), { 
    #print(all_params)
    
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB
    
    gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    
    etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
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
  
    
    # Erreur au carré
    return(c(eq_Snv, eq_Cr, eq_Cs, eq_Csr))
  })
}

out <- ode(
  y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
        Sv = 0,Csv = 0,Csrv =0,Crv = 0,
        Isnv = 0,Isrnv = 0,Irnv = 0,
        Isv = 0,Isrv = 0,Irv = 0),
  times = seq(from=0,to=80000,by=100),
  func = SCISsrV.model_2,
  parms = c(betaC = 0.0152625495139866,
            f = 0.915346433742047,
            as = 3.46118721459732e-06,
            ar = 3.46118721459732e-06,
            betaI = 0, 
            time_until_recovery_without_ATB_s = 50,
            time_until_recovery_without_ATB_r = 50,
            dps = 98,
            dpr = 98,
            prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
            time_until_decolo_by_bystander_ATB = 10,
            # prob_minority_strain_when_colonised = 0.1,
            prob_specific_exposure = 1,
            time_until_decolo_by_specific_ATB = 14,
            # prob_minority_strain_when_infected = 0.1,
            prob_specific_exposure_r = 1,
            time_until_decolo_by_specific_ATB_r = 14,
            sigmaC1 = 1.2*1e-3, 
            sigmaI1 = 1.2*1e-3,
            vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
            vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
) %>% as.data.frame()

out %>% 
  mutate(testC = Csrnv /(Csrnv + Csnv),
         testI = Isrnv /(Isrnv + Isnv))%>%
  mutate(prop_col = (Csnv + Csrnv + Crnv)/(Snv+Csnv+Csrnv+Crnv+Isnv+Isrnv+Irnv),
         prop_col_res = Crnv/(Csrnv + Csnv + Crnv))%>%
  slice_tail(n = 1)

get_equilibrium(c(betaC = 0.0152625495139866,
                  f = 0.915346433742047,
                  as = 3.46118721459732e-06,
                  ar = 3.46118721459732e-06,
                  betaI = 0, 
                  time_until_recovery_without_ATB_s = 50,
                  time_until_recovery_without_ATB_r = 50,
                  dps = 98,
                  dpr = 98,
                  prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
                  time_until_decolo_by_bystander_ATB = 10,
                  # prob_minority_strain_when_colonised = 0.1,
                  prob_specific_exposure = 1,
                  time_until_decolo_by_specific_ATB = 14,
                  # prob_minority_strain_when_infected = 0.1,
                  prob_specific_exposure_r = 1,
                  time_until_decolo_by_specific_ATB_r = 14,
                  sigmaC1 = 1.2*1e-3, 
                  sigmaI1 = 1.2*1e-3,
                  vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
                  vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0))

out %>% 
  mutate(test = Csrnv /(Csrnv + Csnv))%>%
  pivot_longer(!time, names_to="name", values_to = "value" ) %>%
  ggplot(aes(x = time, y = value, color = name))+
  geom_line()+
  facet_wrap(~name, scales = "free_y")+
  theme_bw()


results_2 <- tibble()

for(betaC_i in seq(0.01,0.02,0.001)){
  for(f_i in seq(0.9, 0.95, 0.01)){
    
    out <- ode(
      y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
            Sv = 0,Csv = 0,Csrv =0,Crv = 0,
            Isnv = 0,Isrnv = 0,Irnv = 0,
            Isv = 0,Isrv = 0,Irv = 0),
      times = seq(from=0,to=1000000,by=100),
      func = SCISsrV.model_2,
      parms = c(betaC = betaC_i,
                f = f_i,
                as = 3.46118721459732e-06,
                ar = 3.46118721459732e-06,
                betaI = 0, 
                time_until_recovery_without_ATB_s = 50,
                time_until_recovery_without_ATB_r = 50,
                dps = 98,
                dpr = 98,
                prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
                time_until_decolo_by_bystander_ATB = 10,
                # prob_minority_strain_when_colonised = 0.1,
                prob_specific_exposure = 1,
                time_until_decolo_by_specific_ATB = 14,
                # prob_minority_strain_when_infected = 0.1,
                prob_specific_exposure_r = 1,
                time_until_decolo_by_specific_ATB_r = 14,
                sigmaC1 = 1.2*1e-3, 
                sigmaI1 = 1.2*1e-3,
                vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
                vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
    ) %>% as.data.frame()
    
    results_2 <- bind_rows(results_2, tibble(betaC = betaC_i, 
                                         f = f_i, 
                                         prop_res =out$Crnv[nrow(out)]/(out$Crnv[nrow(out)]+ out$Csnv[nrow(out)] + out$Csrnv[nrow(out)]),
                                         error = abs(out$Crnv[nrow(out)] - out$Crnv[nrow(out)-1])))
    
  }
}

results_2 %>%
  ggplot(aes(x = betaC, y = f, fill = prop_res))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_res, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()



SCISsrV.model_3 <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    # Total population size
    N=Snv+Csnv+Csrnv+Crnv+Isnv+Isrnv+Irnv+Sv+Csv+Csrv+Crv+Isv+Isrv+Irv
    
    ## Transmission coefficients
    Bs = betaC*(Csnv+Csv+Csrnv+Csrv)/N  + betaI*(Isnv+Isv+Isrnv+Isrv)/N
    Br = betaC*f*(Crnv+Crv)/N + betaI*f*(Irnv+Irv)/N
    Bsv = betaC*(1-vftcs)*(Csnv+Csv+Csrnv+Csrv)/N  + betaI*(1-vftis)*(Isnv+Isv+Isrnv+Isrv)/N
    Brv = betaC*(1-vftcr)*f*(Crnv+Crv)/N + betaI*(1-vftir)*f*(Irnv+Irv)/N
    
    
    # Bystander antibiotic exposure on colonised individuals
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    
    # Bystander antibiotic exposure on infected individuals
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    # Specific antibiotic exposure when infected
    if(!("gammaA_s" %in% names(param))){
      gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB
    }
    if(!("gammaA_r" %in% names(param))){
      gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    }
    
    # Recovery rates (when antibiotics aren't taken)
    etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    # Sigma 2 
    sigmaC2 = lambdaA
    sigmaI2 = lambdaA_I + gammaA_s
    
    # Equations of non-vaccinated compartments
    dSnv <- - Bs*Snv - Br*Snv + 1/dps*Csnv+1/dps*Csrnv+1/dpr*Crnv +lambdaA*Csnv + gammaA_r*Irnv + lambdaA_I*Isnv + gammaA_s*Isnv
    dCsnv <- Bs*Snv - 1/dps*Csnv - lambdaA*Csnv  -as*Csnv + etas*Isnv - sigmaC1*Csnv*Br
    dCsrnv <- sigmaC1*Csnv*Br - 1/dps*Csrnv - as*Csrnv - sigmaC2*Csrnv + etas*Isrnv
    dCrnv <- Br*Snv - 1/dpr*Crnv +  sigmaC2*Csrnv - ar*Crnv + etar*Irnv
    dIsnv <- as*Csnv - etas*Isnv - lambdaA_I*Isnv - gammaA_s*Isnv 
    dIsrnv <- as*Csrnv - etas*Isrnv  - sigmaI2*Isrnv
    dIrnv <- ar*Crnv - etar*Irnv -gammaA_r*Irnv + sigmaI2*Isrnv
    
    # Equations of vaccinated compartments
    dSv <- - Bsv*Sv - Brv*Sv + 1/dps*1/(1-vfds)*(Csv+Csrv)+1/dpr*1/(1-vfdr)*Crv + lambdaA*Csv + gammaA_r*Irv + lambdaA_I*Isv + gammaA_s*Isv
    dCsv <- Bsv*Sv - 1/dps*1/(1-vfds)*Csv - lambdaA*Csv - as*(1-vfis)*Csv + etas*1/(1-vfrs)*Isv - sigmaC1*Csv*Brv
    dCsrv <- sigmaC1*Csv*Brv - 1/dps*1/(1-vfds)*Csrv - as*(1-vfis)*Csrv - sigmaC2*Csrv + etas*1/(1-vfrs)*Isrv
    dCrv <- Brv*Sv  - 1/dpr*1/(1-vfdr)*Crv- ar*(1-vfir)*Crv + etar*1/(1-vfrr)*Irv + sigmaC2*Csrv
    dIsv <- as*(1-vfis)*Csv - etas*1/(1-vfrs)*Isv   - lambdaA_I*Isv - gammaA_s*Isv 
    dIsrv <- as*(1-vfis)*Csrv - etas*1/(1-vfrs)*Isrv  - sigmaI2*Isrv
    dIrv <- ar*(1-vfir)*Crv - etar*1/(1-vfrr)*Irv - gammaA_r*Irv + sigmaI2*Isrv
    
    
    res<-c(dSnv, dCsnv,dCsrnv, dCrnv,
           dSv, dCsv,dCsrv, dCrv, 
           dIsnv,dIsrnv, dIrnv, 
           dIsv,dIsrv, dIrv
    )
    
    list(res)
    
  })
  
}

out <- ode(
  y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
        Sv = 0,Csv = 0,Csrv =0,Crv = 0,
        Isnv = 0,Isrnv = 0,Irnv = 0,
        Isv = 0,Isrv = 0,Irv = 0),
  times = seq(from=0,to=1000000,by=100),
  func = SCISsrV.model_3,
  parms = c(betaC = 0.01,#0.0152625495139866,
            f = 0.9, #0.94,#0.915346433742047,
            as = 1e-6,#3.46118721459732e-06,
            ar = 1e-6, #3.46118721459732e-06,
            betaI = 0, 
            time_until_recovery_without_ATB_s = 50,
            time_until_recovery_without_ATB_r = 50,
            dps = 98,
            dpr = 98,
            prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
            time_until_decolo_by_bystander_ATB = 10,
            # prob_minority_strain_when_colonised = 0.1,
            prob_specific_exposure = 1,
            time_until_decolo_by_specific_ATB = 14,
            # prob_minority_strain_when_infected = 0.1,
            prob_specific_exposure_r = 1,
            time_until_decolo_by_specific_ATB_r = 14,
            sigmaC1 = 1, 
            sigmaI1 = 1,
            vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
            vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
) %>% as.data.frame()

out %>% 
  mutate(testC = Csrnv /(Csrnv + Csnv),
         testI = Isrnv /(Isrnv + Isnv))%>%
  mutate(prop_col = (Csnv + Csrnv + Crnv)/(Snv+Csnv+Csrnv+Crnv+Isnv+Isrnv+Irnv),
         prop_col_res = Crnv/(Csrnv + Csnv + Crnv))%>%
  slice_tail(n = 1)

out %>% 
  mutate(test = Csrnv /(Csrnv + Csnv))%>%
  pivot_longer(!time, names_to="name", values_to = "value" ) %>%
  ggplot(aes(x = time, y = value, color = name))+
  geom_line()+
  facet_wrap(~name, scales = "free_y")+
  theme_bw()

out <- ode(
  y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
        Sv = 0,Csv = 0,Csrv =0,Crv = 0,
        Isnv = 0,Isrnv = 0,Irnv = 0,
        Isv = 0,Isrv = 0,Irv = 0),
  times = seq(from=0,to=1000000,by=100),
  func = SCISsrV.model_3,
  parms = c(betaC = 0.01,#0.0152625495139866,
            f = 0.9, #0.94,#0.915346433742047,
            as = 1e-6,#3.46118721459732e-06,
            ar = 1e-6, #3.46118721459732e-06,
            betaI = 0, 
            time_until_recovery_without_ATB_s = 50,
            time_until_recovery_without_ATB_r = 50,
            dps = 98,
            dpr = 98,
            prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
            time_until_decolo_by_bystander_ATB = 10,
            # prob_minority_strain_when_colonised = 0.1,
            prob_specific_exposure = 1,
            time_until_decolo_by_specific_ATB = 14,
            # prob_minority_strain_when_infected = 0.1,
            prob_specific_exposure_r = 1,
            time_until_decolo_by_specific_ATB_r = 14,
            sigmaC1 = 1, 
            sigmaI1 = 1,
            vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
            vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
) %>% as.data.frame()

#### Simulation un an sans vaccin et un an avec vaccin ####

params <- list(betaC = 0.0152625495139866,
            f = 0.94,#0.915346433742047,
            as = 3.46118721459732e-06,
            ar = 3.46118721459732e-06,
            betaI = 0, 
            time_until_recovery_without_ATB_s = 50,
            time_until_recovery_without_ATB_r = 50,
            dps = 98,
            dpr = 98,
            prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
            time_until_decolo_by_bystander_ATB = 10,
            # prob_minority_strain_when_colonised = 0.1,
            prob_specific_exposure = 1,
            time_until_decolo_by_specific_ATB = 14,
            # prob_minority_strain_when_infected = 0.1,
            prob_specific_exposure_r = 1,
            time_until_decolo_by_specific_ATB_r = 14,
            sigmaC1 = 1, 
            sigmaI1 = 1,
            eq_Snv = out$Snv[nrow(out)],
            eq_Csnv = out$Csnv[nrow(out)],
            eq_Csrnv = out$Csrnv[nrow(out)],
            eq_Crnv = out$Crnv[nrow(out)],
            eq_Isnv = out$Isnv[nrow(out)],
            eq_Isrnv = out$Isrnv[nrow(out)],
            eq_Irnv = out$Irnv[nrow(out)]
            )

out <- ode(
  y = c(Snv = params$eq_Snv,Csnv = params$eq_Csnv,Csrnv = params$eq_Csrnv,Crnv = params$eq_Crnv,
      Sv = 0,Csv = 0,Csrv = 0, Crv =0, 
      Isnv = params$eq_Isnv ,Isrnv = params$eq_Isrnv ,Irnv = params$eq_Irnv , 
      Isv = 0,Isrv = 0, Irv =0),
  times = seq(from=0,to=370,by=10),
  func = SCISsrV.model_3,
  parms = c(params, vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
            vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
) %>% as.data.frame()

out %>%
  slice_tail(n=1) %>% 
  select(Snv, Csnv, Csrnv, Crnv, Isnv,  Isrnv, Irnv)

results <- tibble()

for(betaC_i in seq(0.01,0.02,0.001)){
  for(f_i in seq(0.9, 0.95, 0.01)){
    
    out <- ode(
      y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
            Sv = 0,Csv = 0,Csrv =0,Crv = 0,
            Isnv = 0,Isrnv = 0,Irnv = 0,
            Isv = 0,Isrv = 0,Irv = 0),
      times = seq(from=0,to=1000000,by=100),
      func = SCISsrV.model_3,
      parms = c(betaC = betaC_i,#0.0152625495139866,
                f = f_i, #0.94,#0.915346433742047,
                as = 3.46118721459732e-06,
                ar = 3.46118721459732e-06,
                betaI = 0, 
                time_until_recovery_without_ATB_s = 50,
                time_until_recovery_without_ATB_r = 50,
                dps = 98,
                dpr = 98,
                prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
                time_until_decolo_by_bystander_ATB = 10,
                # prob_minority_strain_when_colonised = 0.1,
                prob_specific_exposure = 1,
                time_until_decolo_by_specific_ATB = 14,
                # prob_minority_strain_when_infected = 0.1,
                prob_specific_exposure_r = 1,
                time_until_decolo_by_specific_ATB_r = 14,
                sigmaC1 = 1, 
                sigmaI1 = 1,
                vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
                vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
    ) %>% as.data.frame()
    
    results <- bind_rows(results, tibble(betaC = betaC_i, 
                                         f = f_i, 
                                         prop_res =out$Crnv[nrow(out)]/(out$Crnv[nrow(out)]+ out$Csnv[nrow(out)] + out$Csrnv[nrow(out)]),
                                         error = abs(out$Crnv[nrow(out)] - out$Crnv[nrow(out)-1])))

  }
}

results %>%
  ggplot(aes(x = betaC, y = f, fill = prop_res))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_res, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()

SCISsrV.model_4 <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    # Total population size
    N=Snv+Csnv+Csrnv+Crnv+Crsnv+Isnv+Isrnv+Irnv+Irsnv+Sv+Csv+Csrv+Crv+Crsv+Isv+Isrv+Irv+Irsv
    
    ## Transmission coefficients
    Bs = betaC*(Csnv+Csv+Csrnv+Csrv)/N  + betaI*(Isnv+Isv+Isrnv+Isrv)/N
    Br = betaC*f*(Crnv+Crv+Crsnv+Crsv)/N + betaI*f*(Irnv+Irv+Irsnv+Irsv)/N
    Bsv = betaC*(1-vftcs)*(Csnv+Csv+Csrnv+Csrv)/N  + betaI*(1-vftis)*(Isnv+Isv+Isrnv+Isrv)/N
    Brv = betaC*(1-vftcr)*f*(Crnv+Crv+Crsnv+Crsv)/N + betaI*(1-vftir)*f*(Irnv+Irv+Irsnv+Irsv)/N
    
    
    # Bystander antibiotic exposure on colonised individuals
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    
    # Bystander antibiotic exposure on infected individuals
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    # Specific antibiotic exposure when infected
    if(!("gammaA_s" %in% names(param))){
      gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB
    }
    if(!("gammaA_r" %in% names(param))){
      gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    }
    
    # Recovery rates (when antibiotics aren't taken)
    etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    # Sigma 2 
    sigmaC2 = lambdaA
    sigmaI2 = lambdaA_I + gammaA_s
    
    # Equations of non-vaccinated compartments
    dSnv <- - Bs*Snv - Br*Snv + 1/dps*Csnv + 1/dps*Csrnv + 1/dpr*Crnv + 1/dpr*Crsnv+ lambdaA*Csnv + gammaA_r*Irnv +gammaA_r*Irsnv + lambdaA_I*Isnv + gammaA_s*Isnv
    dCsnv <- Bs*Snv - 1/dps*Csnv - lambdaA*Csnv  -as*Csnv + etas*Isnv - sigmaC1*Csnv*Br
    dCsrnv <- sigmaC1*Csnv*Br - 1/dps*Csrnv - as*Csrnv - sigmaC2*Csrnv + etas*Isrnv
    dCrnv <- Br*Snv - 1/dpr*Crnv +  sigmaC2*Csrnv - ar*Crnv + etar*Irnv + sigmaC2*Crsnv - sigmaC1*Crnv*Bs
    dCrsnv <- sigmaC1*Crnv*Bs - 1/dpr*Crsnv - ar*Crsnv - sigmaC2*Crsnv + etar*Irsnv
    dIsnv <- as*Csnv - etas*Isnv - lambdaA_I*Isnv - gammaA_s*Isnv - sigmaI1*Isnv*Br
    dIsrnv <- as*Csrnv - etas*Isrnv  - sigmaI2*Isrnv + sigmaI1*Isnv*Br
    dIrnv <- ar*Crnv - etar*Irnv -gammaA_r*Irnv + sigmaI2*Isrnv - sigmaI1r*Irnv*Bs + lambdaA_I*Irsnv 
    dIrsnv <- ar*Crsnv - etar*Irsnv -gammaA_r*Irsnv + sigmaI1r*Irnv*Bs - lambdaA_I*Irsnv 
    
    # Equations of vaccinated compartments
    dSv <- - Bsv*Sv - Brv*Sv + 1/dps*1/(1-vfds)*(Csv+Csrv)+1/dpr*1/(1-vfdr)*(Crv+Crsv) + lambdaA*Csv + gammaA_r*(Irv+Irsv) + lambdaA_I*Isv + gammaA_s*Isv
    dCsv <- Bsv*Sv - 1/dps*1/(1-vfds)*Csv - lambdaA*Csv - as*(1-vfis)*Csv + etas*1/(1-vfrs)*Isv - sigmaC1*Csv*Brv
    dCsrv <- sigmaC1*Csv*Brv - 1/dps*1/(1-vfds)*Csrv - as*(1-vfis)*Csrv - sigmaC2*Csrv + etas*1/(1-vfrs)*Isrv
    dCrv <- Brv*Sv  - 1/dpr*1/(1-vfdr)*Crv- ar*(1-vfir)*Crv + etar*1/(1-vfrr)*Irv + sigmaC2*(Csrv+Crsv) - sigmaC1*Crv*Bsv 
    dCrsv <- sigmaC1*Crv*Bsv - 1/dpr*1/(1-vfdr)*Crsv - ar*(1-vfir)*Crsv - sigmaC2*Crsv + etar*1/(1-vfrr)*Irsv
    dIsv <- as*(1-vfis)*Csv - etas*1/(1-vfrs)*Isv   - lambdaA_I*Isv - gammaA_s*Isv - sigmaI1*Isv*Brv
    dIsrv <- as*(1-vfis)*Csrv - etas*1/(1-vfrs)*Isrv  - sigmaI2*Isrv + sigmaI1*Isv*Brv
    dIrv <- ar*(1-vfir)*Crv - etar*1/(1-vfrr)*Irv - gammaA_r*Irv + sigmaI2*Isrv - sigmaI1r*Irv*Bsv + lambdaA_I*Irsv 
    dIrsv <- ar*(1-vfir)*Crsv - etar*1/(1-vfrr)*Irsv -gammaA_r*Irsv + sigmaI1r*Irv*Bsv - lambdaA_I*Irsv 
    
    
    ## Cumulative incidence counters
    
    # Infections
    dinccumIsnv = as*Csnv
    dinccumIsrnv = as*Csrnv
    dinccumIrnv = ar*Crnv
    dinccumIrsnv = ar*Crsnv
    dinccumIsv = as*(1-vfis)*Csv
    dinccumIsrv = as*(1-vfis)*Csrv
    dinccumIrv = ar*(1-vfir)*Crv
    dinccumIrsv = ar*(1-vfir)*Crsv
    
    dinccumSelectionOfResistantForNonVaccinated = sigmaI2*Isrnv
    dinccumSelectionOfResistantForVaccinated = sigmaI2*Isrv
    
    
    res<-c(dSnv, dCsnv,dCsrnv, dCrnv,dCrsnv,
           dSv, dCsv,dCsrv, dCrv, dCrsv, 
           dIsnv,dIsrnv, dIrnv, dIrsnv, 
           dIsv,dIsrv, dIrv,dIrsv
    )
    
    list(res, incidenceIsnv = as*Csnv, incidenceIsrnv = as*Csrnv, incidenceIrnv = ar*Crnv,incidenceIrsnv = ar*Crsnv,
         incidenceIsv = as*(1-vfis)*Csv,incidenceIsrv = as*(1-vfis)*Csrv, incidenceIrv = ar*(1-vfir)*Crv, incidenceIrsv = ar*(1-vfir)*Crsv)
    
  })
  
}

results <- tibble()

for(betaC_i in seq(0.01,0.02,0.001)){
  for(f_i in seq(0.9, 0.95, 0.01)){
    
    out <- ode(
      y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000*0.9,Crsnv = 10000*0.1,
            Sv = 0,Csv = 0,Csrv =0,Crv = 0,Crsv = 0,
            Isnv = 0,Isrnv = 0,Irnv = 0,Irsnv = 0,
            Isv = 0,Isrv = 0,Irv = 0,Irsv = 0),
      times = seq(from=0,to=1000000,by=100),
      func = SCISsrV.model_4,
      parms = c(betaC = betaC_i,#0.0152625495139866,
                f = f_i, #0.94,#0.915346433742047,
                as = 3.46118721459732e-06,
                ar = 3.46118721459732e-06,
                betaI = 0, 
                time_until_recovery_without_ATB_s = 50,
                time_until_recovery_without_ATB_r = 50,
                dps = 98,
                dpr = 98,
                prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
                time_until_decolo_by_bystander_ATB = 10,
                # prob_minority_strain_when_colonised = 0.1,
                prob_specific_exposure = 1,
                time_until_decolo_by_specific_ATB = 14,
                # prob_minority_strain_when_infected = 0.1,
                prob_specific_exposure_r = 1,
                time_until_decolo_by_specific_ATB_r = 14,
                sigmaC1 = 1, 
                sigmaI1 = 0,
                sigmaC1r = 1, 
                sigmaI1r = 0,
                vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
                vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
    ) %>% as.data.frame()
    
    results <- bind_rows(results, tibble(betaC = betaC_i, 
                                         f = f_i, 
                                         prop_res =(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)])/(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)]+ out$Csnv[nrow(out)] + out$Csrnv[nrow(out)]),
                                         error = abs(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)] - out$Crnv[nrow(out)-1] - out$Crsnv[nrow(out)-1]),
                                         prop_col = (out$Csnv[nrow(out)]+out$Csrnv[nrow(out)]+out$Crsnv[nrow(out)]+out$Crnv[nrow(out)])/100000
    ))
  }
}

results %>%
  ggplot(aes(x = betaC, y = f, fill = prop_res))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_res, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()

results %>%
  ggplot(aes(x = betaC, y = f, fill = prop_col))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_col, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()


out %>% 
  select(!starts_with("incidence"))%>%
  mutate(test = Csrnv /(Csrnv + Csnv))%>%
  pivot_longer(!time, names_to="name", values_to = "value" ) %>%
  filter(time<1000)%>%
  ggplot(aes(x = time, y = value, color = name))+
  geom_line()+
  facet_wrap(~name, scales = "free_y")+
  theme_bw()

# growth competition
SCISsrV.model_5 <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    # Total population size
    N=Snv+Csnv+Csrnv+Crnv+Crsnv+Isnv+Isrnv+Irnv+Irsnv+Sv+Csv+Csrv+Crv+Crsv+Isv+Isrv+Irv+Irsv
    
    ## Transmission coefficients
    Bs = betaC*(Csnv+Csv+Csrnv+Csrv)/N  + betaI*(Isnv+Isv+Isrnv+Isrv)/N
    Br = betaC*f*(Crnv+Crv+Crsnv+Crsv)/N + betaI*f*(Irnv+Irv+Irsnv+Irsv)/N
    Bsv = betaC*(1-vftcs)*(Csnv+Csv+Csrnv+Csrv)/N  + betaI*(1-vftis)*(Isnv+Isv+Isrnv+Isrv)/N
    Brv = betaC*(1-vftcr)*f*(Crnv+Crv+Crsnv+Crsv)/N + betaI*(1-vftir)*f*(Irnv+Irv+Irsnv+Irsv)/N
    
    
    # Bystander antibiotic exposure on colonised individuals
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    
    # Bystander antibiotic exposure on infected individuals
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB
    
    # Specific antibiotic exposure when infected
    if(!("gammaA_s" %in% names(param))){
      gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB
    }
    if(!("gammaA_r" %in% names(param))){
      gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    }
    
    # Recovery rates (when antibiotics aren't taken)
    etas = (1-prob_specific_exposure) *(1-prob_bystander_exposure)* 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    # Sigma 2 
    sigmaC2 = lambdaA
    sigmaI2 = lambdaA_I + gammaA_s
    
    # Equations of non-vaccinated compartments
    dSnv <- - Bs*Snv - Br*Snv + 1/dps*Csnv + 1/dps*Csrnv + 1/dpr*Crnv + 1/dpr*Crsnv+ lambdaA*Csnv + gammaA_r*Irnv +gammaA_r*Irsnv + lambdaA_I*Isnv + gammaA_s*Isnv
    dCsnv <- Bs*Snv - 1/dps*Csnv - lambdaA*Csnv  -as*Csnv + etas*Isnv - sigmaC1*Csnv*Br + b0*Csrnv
    dCsrnv <- sigmaC1*Csnv*Br - 1/dps*Csrnv - as*Csrnv - sigmaC2*Csrnv + etas*Isrnv - b0*Csrnv + b*Crsnv
    dCrnv <- Br*Snv - 1/dpr*Crnv +  sigmaC2*Csrnv - ar*Crnv + etar*Irnv + sigmaC2*Crsnv - sigmaC1*Crnv*Bs 
    dCrsnv <- sigmaC1*Crnv*Bs - 1/dpr*Crsnv - ar*Crsnv - sigmaC2*Crsnv + etar*Irsnv - b*Crsnv
    dIsnv <- as*Csnv - etas*Isnv - lambdaA_I*Isnv - gammaA_s*Isnv - sigmaI1*Isnv*Br + b0*Isrnv
    dIsrnv <- as*Csrnv - etas*Isrnv  - sigmaI2*Isrnv + sigmaI1*Isnv*Br - b0*Isrnv + b*Irsnv
    dIrnv <- ar*Crnv - etar*Irnv -gammaA_r*Irnv + sigmaI2*Isrnv - sigmaI1r*Irnv*Bs + lambdaA_I*Irsnv 
    dIrsnv <- ar*Crsnv - etar*Irsnv -gammaA_r*Irsnv + sigmaI1r*Irnv*Bs - lambdaA_I*Irsnv - b*Irsnv
    
    # Equations of vaccinated compartments
    dSv <- - Bsv*Sv - Brv*Sv + 1/dps*1/(1-vfds)*(Csv+Csrv)+1/dpr*1/(1-vfdr)*(Crv+Crsv) + lambdaA*Csv + gammaA_r*(Irv+Irsv) + lambdaA_I*Isv + gammaA_s*Isv
    dCsv <- Bsv*Sv - 1/dps*1/(1-vfds)*Csv - lambdaA*Csv - as*(1-vfis)*Csv + etas*1/(1-vfrs)*Isv - sigmaC1*Csv*Brv + b0*Csrv
    dCsrv <- sigmaC1*Csv*Brv - 1/dps*1/(1-vfds)*Csrv - as*(1-vfis)*Csrv - sigmaC2*Csrv + etas*1/(1-vfrs)*Isrv - b0*Csrv + b*Crsv
    dCrv <- Brv*Sv  - 1/dpr*1/(1-vfdr)*Crv- ar*(1-vfir)*Crv + etar*1/(1-vfrr)*Irv + sigmaC2*(Csrv+Crsv) - sigmaC1*Crv*Bsv 
    dCrsv <- sigmaC1*Crv*Bsv - 1/dpr*1/(1-vfdr)*Crsv - ar*(1-vfir)*Crsv - sigmaC2*Crsv + etar*1/(1-vfrr)*Irsv - b*Crsv
    dIsv <- as*(1-vfis)*Csv - etas*1/(1-vfrs)*Isv   - lambdaA_I*Isv - gammaA_s*Isv - sigmaI1*Isv*Brv + b0*Isrv
    dIsrv <- as*(1-vfis)*Csrv - etas*1/(1-vfrs)*Isrv  - sigmaI2*Isrv + sigmaI1*Isv*Brv - b0*Isrv + b*Irsv
    dIrv <- ar*(1-vfir)*Crv - etar*1/(1-vfrr)*Irv - gammaA_r*Irv + sigmaI2*Isrv - sigmaI1r*Irv*Bsv + lambdaA_I*Irsv 
    dIrsv <- ar*(1-vfir)*Crsv - etar*1/(1-vfrr)*Irsv -gammaA_r*Irsv + sigmaI1r*Irv*Bsv - lambdaA_I*Irsv - b*Irsv
    
    
    ## Cumulative incidence counters
    
    # Infections
    dinccumIsnv = as*Csnv
    dinccumIsrnv = as*Csrnv
    dinccumIrnv = ar*Crnv
    dinccumIrsnv = ar*Crsnv
    dinccumIsv = as*(1-vfis)*Csv
    dinccumIsrv = as*(1-vfis)*Csrv
    dinccumIrv = ar*(1-vfir)*Crv
    dinccumIrsv = ar*(1-vfir)*Crsv
    
    dinccumSelectionOfResistantForNonVaccinated = sigmaI2*Isrnv
    dinccumSelectionOfResistantForVaccinated = sigmaI2*Isrv
    
    
    res<-c(dSnv, dCsnv,dCsrnv, dCrnv,dCrsnv,
           dSv, dCsv,dCsrv, dCrv, dCrsv, 
           dIsnv,dIsrnv, dIrnv, dIrsnv, 
           dIsv,dIsrv, dIrv,dIrsv
    )
    
    list(res, incidenceIsnv = as*Csnv, incidenceIsrnv = as*Csrnv, incidenceIrnv = ar*Crnv,incidenceIrsnv = ar*Crsnv,
         incidenceIsv = as*(1-vfis)*Csv,incidenceIsrv = as*(1-vfis)*Csrv, incidenceIrv = ar*(1-vfir)*Crv, incidenceIrsv = ar*(1-vfir)*Crsv)
    
  })
  
}

results <- tibble()

for(betaC_i in seq(0.01,0.02,0.001)){
  for(b_i in seq(0, 0.01, 0.001)){
    
    out <- ode(
      y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000*0.9,Crsnv = 10000*0.1,
            Sv = 0,Csv = 0,Csrv =0,Crv = 0,Crsv = 0,
            Isnv = 0,Isrnv = 0,Irnv = 0,Irsnv = 0,
            Isv = 0,Isrv = 0,Irv = 0,Irsv = 0),
      times = seq(from=0,to=1000000,by=100),
      func = SCISsrV.model_5,
      parms = c(betaC = betaC_i,#0.0152625495139866,
                f = 1, #0.94,#0.915346433742047,
                as = 3.46118721459732e-06,
                ar = 3.46118721459732e-06,
                betaI = 0, 
                time_until_recovery_without_ATB_s = 50,
                time_until_recovery_without_ATB_r = 50,
                dps = 98,
                dpr = 98,
                prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
                time_until_decolo_by_bystander_ATB = 10,
                # prob_minority_strain_when_colonised = 0.1,
                prob_specific_exposure = 1,
                time_until_decolo_by_specific_ATB = 14,
                # prob_minority_strain_when_infected = 0.1,
                prob_specific_exposure_r = 1,
                time_until_decolo_by_specific_ATB_r = 14,
                sigmaC1 = 1, 
                sigmaI1 = 0,
                sigmaC1r = 1, 
                sigmaI1r = 0,
                b = b_i,
                b0 = 4*b_i,
                vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
                vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
    ) %>% as.data.frame()
    
    results <- bind_rows(results, tibble(betaC = betaC_i, 
                                         b = b_i, 
                                         prop_res =(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)])/(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)]+ out$Csnv[nrow(out)] + out$Csrnv[nrow(out)]),
                                         error = abs(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)] - out$Crnv[nrow(out)-1] - out$Crsnv[nrow(out)-1]),
                                         prop_col = (out$Csnv[nrow(out)]+out$Csrnv[nrow(out)]+out$Crsnv[nrow(out)]+out$Crnv[nrow(out)])/100000
                                         ))
    
  }
}

results %>%
  ggplot(aes(x = betaC, y = b, fill = prop_res))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_res, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()

results %>%
  ggplot(aes(x = betaC, y = b, fill = prop_col))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_col, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()

out <- ode(
  y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000*0.9,Crsnv = 10000*0.1,
        Sv = 0,Csv = 0,Csrv =0,Crv = 0,Crsv = 0,
        Isnv = 0,Isrnv = 0,Irnv = 0,Irsnv = 0,
        Isv = 0,Isrv = 0,Irv = 0,Irsv = 0),
  times = seq(from=0,to=1000000,by=100),
  func = SCISsrV.model_5,
  parms = c(betaC = 0.0152625495139866,
            f = 1, #0.94,#0.915346433742047,
            as = 3.46118721459732e-06,
            ar = 3.46118721459732e-06,
            betaI = 0, 
            time_until_recovery_without_ATB_s = 50,
            time_until_recovery_without_ATB_r = 50,
            dps = 98,
            dpr = 98,
            prob_bystander_exposure = (0.06+4.7)/1000, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
            time_until_decolo_by_bystander_ATB = 10,
            # prob_minority_strain_when_colonised = 0.1,
            prob_specific_exposure = 1,
            time_until_decolo_by_specific_ATB = 14,
            # prob_minority_strain_when_infected = 0.1,
            prob_specific_exposure_r = 1,
            time_until_decolo_by_specific_ATB_r = 14,
            sigmaC1 = 1, 
            sigmaI1 = 0,
            sigmaC1r = 1, 
            sigmaI1r = 0,
            b = 0.00215,
            b0 = 4*0.00215,
            vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
            vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
) %>% as.data.frame()

out %>% 
  select(!starts_with("incidence"))%>%
  mutate(test = Csrnv /(Csrnv + Csnv))%>%
  pivot_longer(!time, names_to="name", values_to = "value" ) %>%
  #filter(time<1000)%>%
  ggplot(aes(x = time, y = value, color = name))+
  geom_line()+
  facet_wrap(~name, scales = "free_y")+
  theme_bw()

out %>% 
  mutate(testC = Csrnv /(Csrnv + Csnv),
         testI = Isrnv /(Isrnv + Isnv))%>%
  mutate(prop_col = (Csnv + Csrnv + Crnv+ Crsnv)/(Snv+Csnv+Csrnv+Crnv+Crsnv+Isnv+Isrnv+Irnv+Irsnv),
         prop_col_res = (Crnv+Crsnv)/(Csrnv + Csnv + Crnv+ Crsnv))%>%
  slice_tail(n = 1)


#### GRaphic tils for E. coli  ####

model = SCISsrV.model_5

results <- tibble()

for(betaC_i in seq(0.06,0.07,0.001)){
  # for(f_i in seq(0.82, 0.87, 0.01)){
  for(b_i in seq(0, 0.001, 0.0001)){
    
    out <- ode(
      y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000*0.9,Crsnv = 10000*0.1,
            Sv = 0,Csv = 0,Csrv =0,Crv = 0,Crsv = 0,
            Isnv = 0,Isrnv = 0,Irnv = 0,Irsnv = 0,
            Isv = 0,Isrv = 0,Irv = 0,Irsv = 0),
      times = seq(from=0,to=1000000,by=100),
      func = model,
      parms = c(betaC = betaC_i,#0.0152625495139866,
                f = f_i, #0.94,#0.915346433742047,
                as = 0.0001684211,
                ar = 0.0001684211,
                betaI = betaC_i, 
                time_until_recovery_without_ATB_s = 28,
                time_until_recovery_without_ATB_r = 28,
                dps = 365,
                dpr = 365,
                prob_bystander_exposure = (0.96+0.77)/1000,
                time_until_decolo_by_bystander_ATB = 10,
                # prob_minority_strain_when_colonised = 0.1,
                prob_specific_exposure = 0.8,
                # prob_minority_strain_when_infected = 0.1,
                prob_specific_exposure_r = 0.8,
                gammaA_s = 0.6/2+0.2/7, 
                # psiA = 0.2/7*0.1, 
                gammaA_r = 0.6/2+0.2/9, 
                sigmaC1 = 1, 
                sigmaI1 = 0,
                sigmaC1r = 1, 
                sigmaI1r = 0,
                vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
                vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0,
                b = b_i,
                b0 = 4*b_i)  
    ) %>% as.data.frame()
    
    results <- bind_rows(results, tibble(betaC = betaC_i, 
                                         # f = f_i, 
                                         b = b_i,
                                         prop_res =(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)])/(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)]+ out$Csnv[nrow(out)] + out$Csrnv[nrow(out)]),
                                         error = abs(out$Crnv[nrow(out)]+out$Crsnv[nrow(out)] - out$Crnv[nrow(out)-1] - out$Crsnv[nrow(out)-1]),
                                         prop_col = (out$Csnv[nrow(out)]+out$Csrnv[nrow(out)]+out$Crsnv[nrow(out)]+out$Crnv[nrow(out)])/100000
    ))
  }
}

results %>%
  ggplot(aes(x = betaC, y = f, fill = prop_res))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_res, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()

results %>%
  ggplot(aes(x = betaC, y = f, fill = prop_col))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_col, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()

results %>%
  ggplot(aes(x = betaC, y = b, fill = prop_res))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_res, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()

results %>%
  ggplot(aes(x = betaC, y = b, fill = prop_col))+
  geom_tile()+
  geom_text(mapping = aes(label = round(prop_col, 2), color = error < 1e-10))+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"))+
  theme_minimal()

