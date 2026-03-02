SCISsrV.model <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    # Total population size
    N=Snv+Csnv+Crnv+Isnv+Irnv+Sv+Csv+Crv+Isv+Irv
    
    ## Transmission coefficients
    Bs = betaC*(Csnv+Csv)/N  + betaI*(Isnv + Isv)/N
    Br = betaC*f*(Crnv+Crv)/N + betaI*f*(Irnv+Irv)/N
    Bsv = betaC*(1-vftcs)*(Csnv+Csv)/N  + betaI*(1-vftis)*(Isnv + Isv)/N
    Brv = betaC*(1-vftcr)*f*(Crnv+Crv)/N + betaI*(1-vftir)*f*(Irnv+Irv)/N
    
    # Bystander antibiotic exposure on colonised individuals
    lambdaA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_colonised)
    phiA = prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_colonised
    
    # Bystander antibiotic exposure on infected individuals
    lambdaA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*(1-prob_minority_strain_when_infected)
    phiA_I = (1-prob_specific_exposure) * prob_bystander_exposure * 1/time_until_decolo_by_bystander_ATB*prob_minority_strain_when_infected
    
    # Specific antibiotic exposure when infected
    if(!("gammaA_s" %in% names(param))){
      gammaA_s = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*(1-prob_minority_strain_when_infected)
    }
    if(!("psiA" %in% names(param))){
      psiA = prob_specific_exposure * 1/time_until_decolo_by_specific_ATB*prob_minority_strain_when_infected
    }
    if(!("gammaA_r" %in% names(param))){
      gammaA_r = prob_specific_exposure_r * 1/time_until_decolo_by_specific_ATB_r
    }
    
    # Recovery rates (when antibiotics aren't taken)
    etas = (1-prob_specific_exposure) * 1/time_until_recovery_without_ATB_s
    etar = (1-prob_specific_exposure_r) * 1/time_until_recovery_without_ATB_r
    
    # Equations of non-vaccinated compartments
    dSnv <- - Bs*Snv - Br*Snv + 1/dps*Csnv+1/dpr*Crnv +lambdaA*Csnv + gammaA_r*Irnv + lambdaA_I*eps*Isnv + gammaA_s*Isnv
    dCsnv <- Bs*Snv - 1/dps*Csnv - lambdaA*Csnv - phiA*Csnv - Br*Csnv*thetasr + Bs*Crnv*thetars -as*Csnv + etas*Isnv + lambdaA_I*(1-eps)*Isnv
    dIsnv <- as*Csnv - etas*Isnv - lambdaA_I*(1-eps)*Isnv - lambdaA_I*eps*Isnv - gammaA_s*Isnv - psiA*Isnv - phiA_I*Isnv
    dCrnv <- Br*Snv - 1/dpr*Crnv + phiA*Csnv + Br*Csnv*thetasr - Bs*Crnv*thetars - ar*Crnv + etar*Irnv
    dIrnv <- ar*Crnv - etar*Irnv + psiA*Isnv + phiA_I*Isnv -gammaA_r*Irnv 
    
    # Equations of vaccinated compartments
    dSv <- - Bsv*Sv - Brv*Sv + 1/dps*1/(1-vfds)*Csv+1/dpr*1/(1-vfdr)*Crv + lambdaA*Csv + gammaA_r*Irv + lambdaA_I*eps*Isv + gammaA_s*Isv 
    dCsv <- Bsv*Sv - 1/dps*1/(1-vfds)*Csv - lambdaA*Csv - phiA*Csv - Brv*Csv*thetasr + Bsv*Crv*thetars - as*(1-vfis)*Csv + etas*1/(1-vfrs)*Isv + lambdaA_I*(1-eps)*Isv
    dIsv <- as*(1-vfis)*Csv - etas*1/(1-vfrs)*Isv - lambdaA_I*(1-eps)*Isv - psiA*Isv -phiA_I*Isv - lambdaA_I*eps*Isv - gammaA_s*Isv 
    dCrv <- Brv*Sv  - 1/dpr*1/(1-vfdr)*Crv  + phiA*Csv + Brv*Csv*thetasr - Bsv*Crv*thetars - ar*(1-vfir)*Crv + etar*1/(1-vfrr)*Irv
    dIrv <- ar*(1-vfir)*Crv - etar*1/(1-vfrr)*Irv + psiA*Isv + phiA_I*Isv - gammaA_r*Irv
    
    ## Cumulative incidence counters
    # Colonisations
    dinccumCsnv = Bs*Snv
    dinccumCrnv = Br*Snv
    dinccumCsv = Bsv*Sv
    dinccumCrv = Brv*Sv 
    
    # Infections
    dinccumIsnv = as*Csnv
    dinccumIrnv = ar*Crnv
    dinccumIsv = as*(1-vfis)*Csv
    dinccumIrv = ar*(1-vfir)*Crv
    
    # Natural recoverys
    dinccumSensitiveNaturalRecoveryForNonVaccinated = etas*Isnv
    dinccumResistantNaturalRecoveryForNonVaccinated = etar*Irnv
    dinccumSensitiveNaturalRecoveryForVaccinated = etas*1/(1-vfrs)*Isv
    dinccumResistantNaturalRecoveryForVaccinated = etar*1/(1-vfrr)*Irv
    
    # Recovery by antibiotic for sensitive strain (specific + bystander)
    dinccumSensitiveSpecificAntibioRecoveryForNonVaccinated = gammaA_s*Isnv
    dinccumSensitiveSpecificAntibioRecoveryForVaccinated = gammaA_s*Isv
    dinccumSensitiveBystanderAntibioRecoveryForNonVaccinated = lambdaA_I*Isnv
    dinccumSensitiveBystanderAntibioRecoveryForVaccinated = lambdaA_I*Isv
    
    # Selection of resistance by antibiotic (specific + bystander)
    dinccumSelectionOfResistantBySpecificAntibioForNonVaccinated = psiA*Isnv
    dinccumSelectionOfResistantBySpecificAntibioForVaccinated = psiA*Isv
    dinccumSelectionOfResistantByBystanderAntibioForNonVaccinated = phiA_I*Isnv
    dinccumSelectionOfResistantByBystanderAntibioForVaccinated = phiA_I*Isv
    
    # Recovery by antibiotic for resistant strain
    dinccumResistantSpecificAntibioRecoveryForNonVaccinated = gammaA_r*Irnv
    dinccumResistantSpecificAntibioRecoveryForVaccinated = gammaA_r*Irv
    
    # Super infections 
    dinccumSuperInfectionByResistantForNonVaccinated = Br*Isnv*thetasr
    dinccumSuperInfectionBySensitiveForNonVaccinated = Bs*Irnv*thetars
    dinccumSuperInfectionByResistantForVaccinated = Brv*Isv*thetasr
    dinccumSuperInfectionBySensitiveForVaccinated = Bsv*Irv*thetars
    
    
    res<-c(dSnv, dCsnv, dCrnv,dSv, dCsv, dCrv, dIsnv, dIrnv, dIsv, dIrv, 
           dinccumCsnv, dinccumCrnv, dinccumCsv, dinccumCrv,
           dinccumIsnv, dinccumIrnv, dinccumIsv, dinccumIrv,
           
           dinccumSensitiveNaturalRecoveryForNonVaccinated,dinccumResistantNaturalRecoveryForNonVaccinated,
           dinccumSensitiveNaturalRecoveryForVaccinated,dinccumResistantNaturalRecoveryForVaccinated,
           
           dinccumSensitiveSpecificAntibioRecoveryForNonVaccinated,dinccumSensitiveSpecificAntibioRecoveryForVaccinated,
           dinccumSensitiveBystanderAntibioRecoveryForNonVaccinated,dinccumSensitiveBystanderAntibioRecoveryForVaccinated,
           
           dinccumSelectionOfResistantBySpecificAntibioForNonVaccinated,dinccumSelectionOfResistantBySpecificAntibioForVaccinated,
           dinccumSelectionOfResistantByBystanderAntibioForNonVaccinated,dinccumSelectionOfResistantByBystanderAntibioForVaccinated,
           
           dinccumResistantSpecificAntibioRecoveryForNonVaccinated,dinccumResistantSpecificAntibioRecoveryForVaccinated,
           
           dinccumSuperInfectionByResistantForNonVaccinated,dinccumSuperInfectionBySensitiveForNonVaccinated,
           dinccumSuperInfectionByResistantForVaccinated,dinccumSuperInfectionBySensitiveForVaccinated
    )
    
    list(res,incidenceCsnv = Bs*Snv, incidenceCrnv = Br*Snv,
         incidenceCsv = Bsv*Sv, incidenceCrv = Brv*Sv, 
         incidenceSensitiveNaturalRecoveryForNonVaccinated = etas*Isnv,
         incidenceResistantNaturalRecoveryForNonVaccinated = etar*Irnv,
         incidenceSensitiveNaturalRecoveryForVaccinated = etas*1/(1-vfrs)*Isv,
         incidenceResistantNaturalRecoveryForVaccinated = etar*1/(1-vfrr)*Irv,
         
         incidenceSensitiveSpecificAntibioRecoveryForNonVaccinated = gammaA_s*Isnv,
         incidenceSensitiveSpecificAntibioRecoveryForVaccinated = gammaA_s*Isv,
         incidenceSensitiveBystanderAntibioRecoveryForNonVaccinated = lambdaA_I*Isnv,
         incidenceSensitiveBystanderAntibioRecoveryForVaccinated = lambdaA_I*Isv,
         
         incidenceSelectionOfResistantBySpecificAntibioForNonVaccinated = psiA*Isnv,
         incidenceSelectionOfResistantBySpecificAntibioForVaccinated = psiA*Isv,
         incidenceSelectionOfResistantByBystanderAntibioForNonVaccinated = phiA_I*Isnv,
         incidenceSelectionOfResistantByBystanderAntibioForVaccinated = phiA_I*Isv,
         
         incidenceResistantSpecificAntibioRecoveryForNonVaccinated = gammaA_r*Irnv,
         incidenceResistantSpecificAntibioRecoveryForVaccinated = gammaA_r*Irv,
         
         incidenceSuperInfectionByResistantForNonVaccinated = Br*Isnv*thetasr,
         incidenceSuperInfectionBySensitiveForNonVaccinated = Bs*Irnv*thetars,
         incidenceSuperInfectionByResistantForVaccinated = Brv*Isv*thetasr,
         incidenceSuperInfectionBySensitiveForVaccinated = Bsv*Irv*thetars,
         incidenceIsnv = as*Csnv, incidenceIrnv = ar*Crnv,incidenceIsv = as*(1-vfis)*Csv, incidenceIrv = ar*(1-vfir)*Crv)
    
  })
  
}