
# ==================================================
# R0 computations
# ==================================================

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
