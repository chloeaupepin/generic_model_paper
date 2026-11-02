library(dplyr)
library(deSolve)
library(tidyr)
library(ggplot2)

source(here::here("scripts", "model.R"))
source(here::here("scripts","utils", "functions_Csr.R"))

#### S. aureus ####
## Choose bacteria 
# Comment or uncomment 4 line codes below and rerun everything for each bacteria 
# S_aureus
bacteria = "S_aureus"
folder_name = "S_aureus"
file_name  = "S_aureus_params_eu_with_J01CR.csv"
transmission_by_infected = FALSE


#### Define bacteria ####
Bacteria_params = read.csv(here::here("files",file_name)) 

#### Run analysis equilibrium ####

out <- ode(
  y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
        Sv = 0,Csv = 0,Csrv =0,Crv = 0,
        Isnv = 0,Isrnv = 0,Irnv = 0,
        Isv = 0,Isrv = 0,Irv = 0,
        inccumIsnv=0, inccumIsrnv=0, inccumIrnv=0, inccumIsv=0,inccumIsrv=0, inccumIrv=0,
        inccumSelectionOfResistantForNonVaccinated=0,inccumSelectionOfResistantForVaccinated=0),
  times = seq(from=0,to=1000000,by=100),
  func = SCISsrV.model_3,
  parms = c(betaC = Bacteria_params$betaC,#0.0152625495139866,
            f = 0.94, #0.94,#0.915346433742047,
            as = Bacteria_params$as,#3.46118721459732e-06,
            ar = Bacteria_params$ar, #3.46118721459732e-06,
            betaI = 0, 
            time_until_recovery_without_ATB_s = Bacteria_params$time_until_recovery_without_ATB_s,
            time_until_recovery_without_ATB_r = Bacteria_params$time_until_recovery_without_ATB_r,
            dps = Bacteria_params$dps,
            dpr = Bacteria_params$dpr,
            prob_bystander_exposure = Bacteria_params$prob_bystander_exposure, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
            time_until_decolo_by_bystander_ATB = Bacteria_params$time_until_decolo_by_bystander_ATB,
            # prob_minority_strain_when_colonised = 0.1,
            prob_specific_exposure = Bacteria_params$prob_specific_exposure,
            time_until_decolo_by_specific_ATB = Bacteria_params$time_until_decolo_by_specific_ATB,
            # prob_minority_strain_when_infected = 0.1,
            prob_specific_exposure_r = Bacteria_params$prob_specific_exposure_r,
            time_until_decolo_by_specific_ATB_r = Bacteria_params$time_until_decolo_by_specific_ATB_r,
            sigmaC1 = 1, 
            sigmaI1 = 1,
            vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
            vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
) %>% as.data.frame()

load(here::here("files","S_aureus","equilibrium_results.RData"))

eq_results_to_plot <- eq_results%>%
  select(eq_Snv, eq_Csnv, eq_Crnv, eq_Isnv, eq_Irnv)%>%
  rename(Snv = eq_Snv,Csnv = eq_Csnv, Crnv = eq_Crnv,Isnv = eq_Isnv, Irnv = eq_Irnv )%>%
  pivot_longer(cols = everything(), names_to = "name",values_to = "value")%>%
  mutate(name = factor(name, levels = c("Snv","Csnv","Crnv","Isnv","Irnv")))%>%
  group_by(name)%>%
  summarise(median = median(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975),
            .groups = "drop")

out %>%
  mutate(Csnv = Csnv + Csrnv, Isnv = Isnv + Isrnv)%>%
  select(time, Snv, Csnv, Crnv, Isnv, Irnv)%>%
  pivot_longer(cols = -time, names_to = "name",values_to = "value")%>%
  mutate(name = factor(name, levels = c("Snv","Csnv","Crnv","Isnv","Irnv")))%>%
  ggplot()+
  geom_line(mapping = aes(x = time, y = value, color = name))+
  geom_hline(data = eq_results_to_plot, mapping = aes(yintercept = median, color = name), linetype = "dashed")+
  geom_rect(data = eq_results_to_plot, mapping = aes(xmin = 0, xmax = 1000000,ymin = q025, ymax = q975, fill = name),alpha = 0.2) +
  facet_wrap(~name, scales="free_y", nrow = 1)+
  labs(x = "Time",y = "Number")+
  theme_bw()+theme(legend.position = "none")+
  ylim(0,NA)

ggsave("eq_S_aureus.png", width = 12, height = 3 )

#### E. coli ####
# E_coli
bacteria = "E_coli"
folder_name = "E_coli"
file_name = "E_coli_params.csv"
transmission_by_infected = TRUE # if transmission is allowed, it equals transmission by colonized individuals

#### Define bacteria ####
Bacteria_params = read.csv(here::here("files",file_name)) 

#### Run analysis equilibrium ####

out <- ode(
  y = c(Snv = 80000,Csnv = 10000*0.9,Csrnv = 10000*0.1,Crnv = 10000,
        Sv = 0,Csv = 0,Csrv =0,Crv = 0,
        Isnv = 0,Isrnv = 0,Irnv = 0,
        Isv = 0,Isrv = 0,Irv = 0,
        inccumIsnv=0, inccumIsrnv=0, inccumIrnv=0, inccumIsv=0,inccumIsrv=0, inccumIrv=0,
        inccumSelectionOfResistantForNonVaccinated=0,inccumSelectionOfResistantForVaccinated=0),
  times = seq(from=0,to=1000000,by=100),
  func = SCISsrV.model_3,
  parms = c(betaC = Bacteria_params$betaC,#0.0152625495139866,
            f = Bacteria_params$f, #0.94,#0.915346433742047,
            as = Bacteria_params$as,#3.46118721459732e-06,
            ar = Bacteria_params$ar, #3.46118721459732e-06,
            betaI = 0, 
            time_until_recovery_without_ATB_s = Bacteria_params$time_until_recovery_without_ATB_s,
            time_until_recovery_without_ATB_r = Bacteria_params$time_until_recovery_without_ATB_r,
            dps = Bacteria_params$dps,
            dpr = Bacteria_params$dpr,
            prob_bystander_exposure = Bacteria_params$prob_bystander_exposure, #0.06/1000 #value 0.06 only for S_aureus_params_eu.csv
            time_until_decolo_by_bystander_ATB = Bacteria_params$time_until_decolo_by_bystander_ATB,
            # prob_minority_strain_when_colonised = 0.1,
            prob_specific_exposure = Bacteria_params$prob_specific_exposure,
            #time_until_decolo_by_specific_ATB = Bacteria_params$time_until_decolo_by_specific_ATB,
            # prob_minority_strain_when_infected = 0.1,
            prob_specific_exposure_r = Bacteria_params$prob_specific_exposure_r,
            time_until_decolo_by_specific_ATB_r = Bacteria_params$time_until_decolo_by_specific_ATB_r,
            sigmaC1 = 1.1e-1, 
            sigmaI1 = 1.1e-1,
            gammaA_s = Bacteria_params$gammaA_s,
            psiA = Bacteria_params$psiA,
            gammaA_r = Bacteria_params$gammaA_r,
            vftcs = 0, vftis = 0, vfds = 0, vfis = 0, vfrs = 0,
            vftcr = 0, vftir = 0, vfdr = 0, vfir = 0, vfrr = 0)  
) %>% as.data.frame()

load(here::here("files","E_coli","equilibrium_results.RData"))

eq_results_to_plot <- eq_results%>%
  select(eq_Snv, eq_Csnv, eq_Crnv, eq_Isnv, eq_Irnv)%>%
  rename(Snv = eq_Snv,Csnv = eq_Csnv, Crnv = eq_Crnv,Isnv = eq_Isnv, Irnv = eq_Irnv )%>%
  pivot_longer(cols = everything(), names_to = "name",values_to = "value")%>%
  mutate(name = factor(name, levels = c("Snv","Csnv","Crnv","Isnv","Irnv")))%>%
  group_by(name)%>%
  summarise(median = median(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975),
            .groups = "drop")

out %>%
  mutate(Csnv = Csnv + Csrnv, Isnv = Isnv + Isrnv)%>%
  select(time, Snv, Csnv, Crnv, Isnv, Irnv)%>%
  pivot_longer(cols = -time, names_to = "name",values_to = "value")%>%
  mutate(name = factor(name, levels = c("Snv","Csnv","Crnv","Isnv","Irnv")))%>%
  ggplot()+
  geom_line(mapping = aes(x = time, y = value, color = name))+
  geom_hline(data = eq_results_to_plot, mapping = aes(yintercept = median, color = name), linetype = "dashed")+
  geom_rect(data = eq_results_to_plot, mapping = aes(xmin = 0, xmax = 1000000,ymin = q025, ymax = q975, fill = name),alpha = 0.2) +
  facet_wrap(~name, scales="free_y", nrow=1)+
  labs(x = "Time",y = "Number")+
  theme_bw()+theme(legend.position = "none")+
  ylim(0,NA)

ggsave("eq_E_coli.png", width = 12, height = 3)

