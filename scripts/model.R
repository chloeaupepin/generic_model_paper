SCISsrV.model <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    #vftcr = vftcs
    #vfir = vfis
    #vfdr = vfds
    #vfrr = vfrs
    #dpr = dps
    
    N=Snv+Csnv+Crnv+Isnv+Irnv+Sv+Csv+Crv+Isv+Irv
    #dir = Inf
    
    Bs = betaC*(Csnv+Csv)/N  + betaI*(Isnv + Isv)/N
    Br = betaC*f*(Crnv+Crv)/N + betaI*f*(Irnv+Irv)/N
    Bsv = betaC*(1-vftcs)*(Csnv+Csv)/N  + betaI*(1-vftis)*(Isnv + Isv)/N
    Brv = betaC*(1-vftcr)*f*(Crnv+Crv)/N + betaI*(1-vftir)*f*(Irnv+Irv)/N
    
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
    
    #if(t==100){
    #U = N / betaC * (1/dps + lambdaA + phiA + as *(lambdaA+gammaA_s + psiA)/(etas + lambdaA + gammaA_s + psiA))
    #print(U <= N)
    #E = (betaC * f * U / N - 1/dps - as + as*etas/(etas+gammaA_r))/(etas * psiA * as /((etas + gammaA_r)*(etas + lambdaA + gammaA_s + psiA)) + phiA)
    #print(E <= 0)
    # print(1 + as/(etas + gammaA_r) - E + as/(etas + lambdaA + gammaA_s + psiA) * (-1 + psiA / (etas + gammaA_r) *E) >=0)
    # 
    # print("valeurs equilibres betaI = 0, pas de vaccination, avec suracquisition pour la colonisation, pas d'antibio")
    # print(thetars * (1-f)/f - 1/f*(1+as/(etas+gammaA_r))+1 + as/etas)
    # Cs_star = (N - N / (betaC*f) * (1/dps + as + etas * as /(etas + gammaA_r))-(1+as/(etas+gammaA_r))*N/(betaC*thetars*(f-1))*(1/f * (1/dps + as/(dir*etas + 1))-1/dps))/(thetars * (1-f)/f - 1/f*(1+as/(etas+gammaA_r))+1 + as/etas)
    # Cr_star = N / (betaC*thetars*(f-1)) * (1/f * (1/dps + as/(dir*etas+1) + betaC*thetars*Cs_star/N*(1-f))-1/dps)
    # print(Cs_star)
    # print(Cr_star)
    # print( 1/dps >= Cr_star*thetars*betaC/N*(1-f))
    # print(1/dps >= 1/f * (1/dps + as/(dir*etas+1) + betaC*thetars*Cs_star/N*(1-f)))
    # 
    # print((N-N/(betaC*dps))/(1+as/etas))
    
    # U_star = N/betaC*(1/dps + lambdaA + as + phiA - etas*as/(etas+lambdaA_I+gammaA_s+psiA))
    # E = (betaC * f * U_star / N - 1/dpr - ar + ar*etar/(etar+gammaA_r))/(etar * psiA * as /((etar + gammaA_r)*(etas + lambdaA_I + gammaA_s + psiA)) + phiA)
    # G = 1+ar/(etar+gammaA_r)-E*(1+as/(etas+lambdaA_I+gammaA_s+psiA)*(1+psiA/(etar+gammaA_r)))
    # Cr_star = (N-U_star)/G
    # Cs_star = - E * Cr_star
    # Is_star = as/(etas+lambdaA_I+gammaA_s+psiA)*Cs_star
    # Ir_star = (ar*Cr_star+psiA*Is_star)/(etar+gammaA_r)
    # 
    # J = matrix(c(betaC/N*(N-2*Cs_star-Cr_star-Is_star-Ir_star)-1/dps-lambdaA-as-phiA, -betaC/N*Cs_star,-betaC/N*Cs_star+etas, -betaC/N*Cs_star,
    #              -betaC/N*f*Cr_star+phiA, betaC*f/N*(N-Cs_star-2*Cr_star-Is_star-Ir_star)-1/dpr-ar, -betaC*f/N*Cr_star, - betaC/N*f*Cr_star+etar,
    #              as, 0, -etas-lambdaA_I-gammaA_s-psiA,0,
    #              0,ar,psiA,- etar - gammaA_r), nrow=4, byrow = TRUE)
    # 
    # print(eigen(J))
    # U = N / betaC * (1/dps + lambdaA + phiA + as - etas*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))
    # print(U <= N)
    # E = (betaC * f * U / N - 1/dpr - ar + ar*etar/(etar+gammaA_r))/(etar * (psiA + phiA_I) * as/((etar + gammaA_r)*(etas + lambdaA_I + gammaA_s + psiA + phiA_I)) + phiA)
    # print(E <= 0)
    # 
    # v_star_s = (betaC - 1/dps - lambdaA - phiA - as + etas*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + 1/dps - as + etas*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))
    # print(v_star_s)
    # v_star_r = (betaC*f - 1/dpr - ar + etar*ar/(etar + gammaA_r))/(betaC*f + 1/dpr - ar + etar*ar/(etar + gammaA_r))
    # print(v_star_r)
    # 
    # print(betaC*(1-vftcs) -1/dps*1/(1-vfds) - lambdaA - as*(1-vfis)+etas*1/(1-vfrs)*as*(1-vfis)/(etas*1/(1-vfrs)+lambdaA_I+gammaA_s+psiA+phiA_I) - phiA)
    # print(betaC*f*(1-vftcr) - 1/dpr*1/(1-vfdr)-ar*(1-vfir)+etar*1/(1-vfrr)*ar*(1-vfir)/(etar*1/(1-vfrr)+gammaA_r))
    # 
    # U_vftc = U/(1-vftcs)
    # E_vftc = (betaC * f * U_vftc * (1-vftcs) / N - 1/dpr - ar + ar*etar/(etar+gammaA_r))/(etar * (psiA + phiA_I) * as/((etar + gammaA_r)*(etas + lambdaA_I + gammaA_s + psiA + phiA_I)) + phiA)
    # G_vftc = 1+ar/(etar+gammaA_r)-E_vftc*(1+as/(etas+lambdaA_I+gammaA_s+psiA)*(1+psiA/(etar+gammaA_r)))
    # Cr_vftc_star = (N-U_vftc)/G_vftc
    # 
    # print(Cr_vftc_star)
    #}
    
    dSnv <- - Bs*Snv - Br*Snv + 1/dps*Csnv+1/dpr*Crnv +lambdaA*Csnv + gammaA_r*Irnv + lambdaA_I*eps*Isnv + gammaA_s*Isnv
    dCsnv <- Bs*Snv - 1/dps*Csnv - lambdaA*Csnv - phiA*Csnv - Br*Csnv*thetasr + Bs*Crnv*thetars -as*Csnv + etas*Isnv + lambdaA_I*(1-eps)*Isnv
    dIsnv <- as*Csnv - etas*Isnv - lambdaA_I*(1-eps)*Isnv - lambdaA_I*eps*Isnv - gammaA_s*Isnv - psiA*Isnv - phiA_I*Isnv
    dCrnv <- Br*Snv - 1/dpr*Crnv + phiA*Csnv + Br*Csnv*thetasr - Bs*Crnv*thetars - ar*Crnv + etar*Irnv
    dIrnv <- ar*Crnv - etar*Irnv + psiA*Isnv + phiA_I*Isnv -gammaA_r*Irnv 
    
    dSv <- - Bsv*Sv - Brv*Sv + 1/dps*1/(1-vfds)*Csv+1/dpr*1/(1-vfdr)*Crv + lambdaA*Csv + gammaA_r*Irv + lambdaA_I*eps*Isv + gammaA_s*Isv 
    dCsv <- Bsv*Sv - 1/dps*1/(1-vfds)*Csv - lambdaA*Csv - phiA*Csv - Brv*Csv*thetasr + Bsv*Crv*thetars - as*(1-vfis)*Csv + etas*1/(1-vfrs)*Isv + lambdaA_I*(1-eps)*Isv
    dIsv <- as*(1-vfis)*Csv - etas*1/(1-vfrs)*Isv - lambdaA_I*(1-eps)*Isv - psiA*Isv -phiA_I*Isv - lambdaA_I*eps*Isv - gammaA_s*Isv 
    dCrv <- Brv*Sv  - 1/dpr*1/(1-vfdr)*Crv  + phiA*Csv + Brv*Csv*thetasr - Bsv*Crv*thetars - ar*(1-vfir)*Crv + etar*1/(1-vfrr)*Irv
    dIrv <- ar*(1-vfir)*Crv - etar*1/(1-vfrr)*Irv + psiA*Isv + phiA_I*Isv - gammaA_r*Irv
    
    dinccumCsnv = Bs*Snv
    dinccumCrnv = Br*Snv
    dinccumCsv = Bsv*Sv
    dinccumCrv = Brv*Sv 
    
    dinccumIsnv = as*Csnv
    dinccumIrnv = ar*Crnv
    dinccumIsv = as*(1-vfis)*Csv
    dinccumIrv = ar*(1-vfir)*Crv
    
    dinccumSensitiveNaturalRecoveryForNonVaccinated = etas*Isnv
    dinccumResistantNaturalRecoveryForNonVaccinated = etar*Irnv
    dinccumSensitiveNaturalRecoveryForVaccinated = etas*1/(1-vfrs)*Isv
    dinccumResistantNaturalRecoveryForVaccinated = etar*1/(1-vfrr)*Irv
    
    dinccumSensitiveSpecificAntibioRecoveryForNonVaccinated = gammaA_s*Isnv
    dinccumSensitiveSpecificAntibioRecoveryForVaccinated = gammaA_s*Isv
    dinccumSensitiveBystanderAntibioRecoveryForNonVaccinated = lambdaA_I*Isnv
    dinccumSensitiveBystanderAntibioRecoveryForVaccinated = lambdaA_I*Isv
    
    dinccumSelectionOfResistantBySpecificAntibioForNonVaccinated = psiA*Isnv
    dinccumSelectionOfResistantBySpecificAntibioForVaccinated = psiA*Isv
    dinccumSelectionOfResistantByBystanderAntibioForNonVaccinated = phiA_I*Isnv
    dinccumSelectionOfResistantByBystanderAntibioForVaccinated = phiA_I*Isv
    
    dinccumResistantSpecificAntibioRecoveryForNonVaccinated = gammaA_r*Irnv
    dinccumResistantSpecificAntibioRecoveryForVaccinated = gammaA_r*Irv
    
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