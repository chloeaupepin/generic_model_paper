library(purrr)
library(tidyr)
library(dplyr)
library(deSolve)
library(furrr)
library(stringr)

source(here::here("scripts", "model.R"))
source(here::here("scripts","utils", "functions.R"))

n_cores = parallel::detectCores() - 1
plan(multisession, workers = n_cores)


set.seed(123)  # global seed for reproducibility

#### Define settings ####

## Population size
population_size = 100000

## Choose number of bacteria to generate
n <- 50

## Choose bacteria 
# S_aureus
bacteria = "S_aureus"
folder_name = "S_aureus"
file_name  = "S_aureus_params10.csv"

Bacteria_params_S_aureus = read.csv(here::here("files",file_name)) 

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

df_bacteria_S_aureus <- add_variability(Bacteria_params_S_aureus, cols_to_noise, n ,0.1)%>%
  mutate(bacteria_name = bacteria)%>%
  mutate(betaI = 0)

df_bacteria_S_aureus <- compute_antibiotic_associated_flux(df_bacteria_S_aureus, T)


# E_coli
bacteria = "E_coli"
folder_name = "E_coli"
file_name = "E_coli_params_primavera4.csv"

Bacteria_params_E_coli = read.csv(here::here("files",file_name)) 

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

df_bacteria_E_coli <- add_variability(Bacteria_params_E_coli, cols_to_noise, n ,0.1)%>%
  mutate(bacteria_name = bacteria) %>%
  mutate(betaI = betaC)

df_bacteria_E_coli <- compute_antibiotic_associated_flux(df_bacteria_E_coli, F)


#### Define bacteria ####

df_bacteria <- bind_rows(df_bacteria_S_aureus, df_bacteria_E_coli)

df_bacteria <- set_resistant_and_sensitive_param_equal(df_bacteria, c("dpr", "ar", "time_until_recovery_without_ATB_r"))


df_bacteria <- df_bacteria %>%
  mutate(prob_specific_exposure_r = pmin(1,prob_specific_exposure_r),
         prob_specific_exposure = pmin(1,prob_specific_exposure),
         prob_minority_strain_when_colonised = pmin(1,prob_minority_strain_when_colonised),
         prob_minority_strain_when_infected = pmin(1,prob_minority_strain_when_infected),
         prob_bystander_exposure = pmin(1,prob_bystander_exposure)
  )

df_bacteria <- df_bacteria %>%
  mutate(g = betaC + betaI*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I) - ( 1/dps + lambdaA + as - etas * as / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA))%>%
  mutate(h = -f + (1/dpr + ar - etar * ar / (etar + gammaA_r)) /(1/dps + lambdaA + as - etas * as / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar/(etar + gammaA_r)))


#### Run analysis equilibrium ####

# Compute equilibrium with analytical expressions 
eq_results <- compute_equilibrium(df_bacteria, population_size)

# Check stability of equilibrium
eq_results <- eq_results %>%
  mutate(stability = pmap(., check_eq_stability)) %>%
  unnest_wider(stability) 


#### Simulation de un an sans vaccin et un an avec vaccin ####
# Simulate 1 year without vaccine for each parameter set at equilibrium

# Choose vaccination scenario
vaccine_scenarios_complete_df <- bind_rows(
  expand.grid("Vperc" = c(0.9),"vftcs" = c(0,0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=0,"vfdr"=0,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0,name = "vftc"),
  expand.grid("Vperc" = c(0.9),"vfds" = c(0,0.3,0.6,0.9)) %>% mutate("vfdr" = vfds, "vftcs"=0, "vftcr"=0,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0,name = "vfd"),
  expand.grid("Vperc" = c(0.9),"vfis" = c(0,0.3,0.6,0.9)) %>% mutate("vfir" = vfis,"vftcs"=0, "vftcr"=0, "vfds"=0, "vfdr"=0, "vfrs"=0, "vfrr"=0, name = "vfi"),
  expand.grid("Vperc" = c(0.9),"vftcs" = c(0,0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfd"),
  expand.grid("Vperc" = c(0.9),"vftcs" = c(0,0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfds"=0,"vfdr"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfi"),
  expand.grid("Vperc" = c(0.9),"vfds" = c(0,0.3,0.6,0.9)) %>% mutate("vfdr" = vfds, "vfis"=vfds,"vfir"=vfds,"vftcs"=0, "vftcr"=0,"vfrs"=0, "vfrr"=0, name = "vfd_vfi"),
  expand.grid("Vperc" = c(0.9),"vftcs" = c(0,0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfrs"=0, "vfrr"=0, name = "vftc_vfd_vfi")
) %>% mutate(vaccine_id = row_number()) 



# Add vaccine parameters to each parameter set
df_for_vaccine_simulations <- add_vaccine_parameters(eq_results, vaccine_scenarios_complete_df)%>%
  mutate(vftis = case_when(bacteria_name == "E_coli" ~ vftcs, T~ 0), vftir = case_when(bacteria_name == "E_coli" ~ vftcr, T~ 0), )

df_for_vaccine_simulations_test <- df_for_vaccine_simulations %>%
  mutate(hv = -f + (1/(dpr*(1-vfdr)) + ar*(1-vfir) - etar * ar*(1-vfir) / (etar + gammaA_r)) /(1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas * as*(1-vfis) / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as*(1-vfis)/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar*(1-vfir)/(etar + gammaA_r))) %>%
  mutate(gv = (1-vftcs)*(betaC + betaI*as*(1-vfis)/(etas + lambdaA_I + gammaA_s + psiA + phiA_I)) - ( 1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas * as*(1-vfis) / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA))%>%
  mutate(
    eq_Sv = N/(betaC*(1-vftcs) + betaI*(1-vftcs)*as*(1-vfis)/His)*(1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas*as*(1-vfis)/His + phiA),
    E1v = f*eq_Sv/N*(betaC*(1-vftcr) + betaI*(1-vftcr)*ar*(1-vfir)/(etar+gammaA_r))-1/(dpr*(1-vfdr))- ar*(1-vfir) + etar*ar*(1-vfir)/(etar+gammaA_r),
    E2v = (betaI*(1-vftcr)*f*eq_Sv/N +etar)*(psiA + phiA_I)/(etar+gammaA_r)*as*(1-vfis)/His +phiA,
    Ev = E1v/E2v,
    Gv = 1+ ar*(1-vfir)/(etar+gammaA_r)-Ev*(1+as*(1-vfis)/His*(1+(psiA + phiA_I)/(etar+gammaA_r))),
    eq_Crv = (N-eq_Sv)/Gv,
    eq_Csv = -Ev*eq_Crv
  )%>%
  mutate(test1 = 1/(1-E),
         test2 = case_when( hv> 0 & gv > 0 ~ 1/(1-Ev), hv > 0 & gv< 0 ~ NA)) %>%
  mutate(test = test2 - test1)

# Simulate one year with vaccine for each parameter set
results <- df_for_vaccine_simulations %>%
  apply_function_on_df(simulate_1y_with_vaccine, "res_1y_wv") %>%
  unnest_wider(res_1y_wv, names_sep = "_")


#### Results manipulation ####
# Compute resistance change
results_with_outputs <- results %>% 
  mutate(hv = -f + (1/(dpr*(1-vfdr)) + ar*(1-vfir) - etar * ar*(1-vfir) / (etar + gammaA_r)) /(1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas * as*(1-vfis) / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as*(1-vfis)/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar*(1-vfir)/(etar + gammaA_r))) %>%
  mutate(gv = (1-vftcs)*(betaC + betaI*as*(1-vfis)/(etas + lambdaA_I + gammaA_s + psiA + phiA_I)) - ( 1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas * as*(1-vfis) / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA))%>%
  mutate(
    eq_Sv = N/(betaC*(1-vftcs) + betaI*(1-vftcs)*as*(1-vfis)/His)*(1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas*as*(1-vfis)/His + phiA),
    E1v = f*eq_Sv/N*(betaC*(1-vftcr) + betaI*(1-vftcr)*ar*(1-vfir)/(etar+gammaA_r))-1/(dpr*(1-vfdr))- ar*(1-vfir) + etar*ar*(1-vfir)/(etar+gammaA_r),
    E2v = (betaI*(1-vftcr)*f*eq_Sv/N +etar)*(psiA + phiA_I)/(etar+gammaA_r)*as*(1-vfis)/His +phiA,
    Ev = E1v/E2v,
    Gv = 1+ ar*(1-vfir)/(etar+gammaA_r)-Ev*(1+as*(1-vfis)/His*(1+(psiA + phiA_I)/(etar+gammaA_r))),
    eq_Crv = (N-eq_Sv)/Gv,
    eq_Csv = -Ev*eq_Crv
  )%>%
  mutate(test1 = 1/(1-E),
         test2 = case_when( hv> 0 & gv > 0 ~ 1/(1-Ev), hv > 0 & gv< 0 ~ NA)) %>%
  mutate(resistance_carriage = res_1y_wv_Crv/population_size, 
         resistance_frequency = res_1y_wv_Crv/(res_1y_wv_Crv + res_1y_wv_Csv))




# Rename vaccine names
results_with_outputs <- results_with_outputs %>%
  mutate(name_renamed= case_when(
    name == "vftc" ~ "va",
    name == "vfd" ~ "vd",
    name == "vfi" ~ "vi",
    name == "vftc_vfd" ~ "va_vd",
    name == "vftc_vfi" ~ "va_vi",
    name == "vfd_vfi" ~ "vd_vi",
    name == "vftc_vfd_vfi" ~ "va_vd_vi"
  )) %>%
  mutate(name_renamed = factor(name_renamed, levels = c("va","vd","vi","va_vd","va_vi","vd_vi","va_vd_vi")))%>%
  mutate(varying_param = case_when(
    name_renamed %in% c("va", "va_vd","va_vi", "va_vd_vi") ~ vftcs,
    name_renamed %in% c("vd", "vd_vi")~ vfds,
    name_renamed %in% c("vi") ~ vfis
  )) %>% 
  mutate(varying_param = factor(varying_param))

# Compute statistics
results_to_plot <- results_with_outputs %>%
  pivot_longer(cols = all_of(c("gv", "hv", "resistance_carriage", "resistance_frequency", "test1","test2")), 
               names_to = "metric_name", values_to = "metric_value")%>%
  group_by(name_renamed, varying_param, metric_name, bacteria_name) %>%
  summarise(median = median(metric_value, na.rm = T),
            q025 = quantile(metric_value, 0.025, na.rm = T),
            q975 = quantile(metric_value, 0.975, na.rm = T),
            .groups = "drop")



library(ggplot2)

ggplot(results_to_plot, aes(x = varying_param, y = median, color = bacteria_name))+
  geom_point()+
  geom_errorbar(aes(ymin = q025, ymax = q975))+
  facet_grid(metric_name ~ name_renamed, scales = "free_y")+
  theme_bw()

ggplot(results_to_plot%>%filter(metric_name %in% c("resistance_frequency", "test1","test2")), aes(x = varying_param, y = median, color = bacteria_name))+
  geom_point()+
  geom_errorbar(aes(ymin = q025, ymax = q975))+
  facet_grid(metric_name ~ name_renamed)+
  theme_bw()


#### au cours du temps ####

simulate_1y_with_vaccine_test <- function(params) {
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
  
  as_tibble(out %>%
    select(time, Snv, Csnv, Crnv, Isnv, Irnv, Sv, Csv, Crv, Isv, Irv,
           inccumIsnv, inccumIrnv, inccumIsv, inccumIrv,
           inccumSelectionOfResistantBySpecificAntibioForNonVaccinated,
           inccumSelectionOfResistantByBystanderAntibioForNonVaccinated,
           inccumSelectionOfResistantBySpecificAntibioForVaccinated,
           inccumSelectionOfResistantByBystanderAntibioForVaccinated))
}

results_test <- df_for_vaccine_simulations %>%
  apply_function_on_df(simulate_1y_with_vaccine_test, "res_1y_wv") %>%
  unnest(res_1y_wv, names_sep = "_")




results_with_outputs <- results_test %>% 
  mutate(hv = -f + (1/(dpr*(1-vfdr)) + ar*(1-vfir) - etar * ar*(1-vfir) / (etar + gammaA_r)) /(1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas * as*(1-vfis) / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as*(1-vfis)/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar*(1-vfir)/(etar + gammaA_r))) %>%
  mutate(gv = (1-vftcs)*(betaC + betaI*as*(1-vfis)/(etas + lambdaA_I + gammaA_s + psiA + phiA_I)) - ( 1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas * as*(1-vfis) / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA))%>%
  mutate(
    eq_Sv = N/(betaC*(1-vftcs) + betaI*(1-vftcs)*as*(1-vfis)/His)*(1/(dps*(1-vfds)) + lambdaA + as*(1-vfis) - etas*as*(1-vfis)/His + phiA),
    E1v = f*eq_Sv/N*(betaC*(1-vftcr) + betaI*(1-vftcr)*ar*(1-vfir)/(etar+gammaA_r))-1/(dpr*(1-vfdr))- ar*(1-vfir) + etar*ar*(1-vfir)/(etar+gammaA_r),
    E2v = (betaI*(1-vftcr)*f*eq_Sv/N +etar)*(psiA + phiA_I)/(etar+gammaA_r)*as*(1-vfis)/His +phiA,
    Ev = E1v/E2v,
    Gv = 1+ ar*(1-vfir)/(etar+gammaA_r)-Ev*(1+as*(1-vfis)/His*(1+(psiA + phiA_I)/(etar+gammaA_r))),
    eq_Crv = (N-eq_Sv)/Gv,
    eq_Csv = -Ev*eq_Crv
  )%>%
  mutate(test1 = 1/(1-E),
         test2 = case_when( hv> 0 & gv > 0 ~ 1/(1-Ev), hv > 0 & gv< 0 ~ NA)) %>%
  mutate(resistance_carriage = res_1y_wv_Crv/population_size, 
         resistance_frequency = res_1y_wv_Crv/(res_1y_wv_Crv + res_1y_wv_Csv),
         resistant_infections = res_1y_wv_inccumIrv/population_size,
         resistant_infections_frequency = res_1y_wv_inccumIrv/(res_1y_wv_inccumIrv+res_1y_wv_inccumIsv))




# Rename vaccine names
results_with_outputs <- results_with_outputs %>%
  mutate(name_renamed= case_when(
    name == "vftc" ~ "va",
    name == "vfd" ~ "vd",
    name == "vfi" ~ "vi",
    name == "vftc_vfd" ~ "va_vd",
    name == "vftc_vfi" ~ "va_vi",
    name == "vfd_vfi" ~ "vd_vi",
    name == "vftc_vfd_vfi" ~ "va_vd_vi"
  )) %>%
  mutate(name_renamed = factor(name_renamed, levels = c("va","vd","vi","va_vd","va_vi","vd_vi","va_vd_vi")))%>%
  mutate(varying_param = case_when(
    name_renamed %in% c("va", "va_vd","va_vi", "va_vd_vi") ~ vftcs,
    name_renamed %in% c("vd", "vd_vi")~ vfds,
    name_renamed %in% c("vi") ~ vfis
  )) %>% 
  mutate(varying_param = factor(varying_param))


# Compute statistics
results_to_plot <- results_with_outputs %>%
  pivot_longer(cols = all_of(c("resistance_carriage", "resistance_frequency", "resistant_infections", "resistant_infections_frequency")), 
               names_to = "metric_name", values_to = "metric_value")%>%
  group_by(res_1y_wv_time, name_renamed, varying_param, metric_name, bacteria_name) %>%
  summarise(median = median(metric_value, na.rm = T),
            q025 = quantile(metric_value, 0.025, na.rm = T),
            q975 = quantile(metric_value, 0.975, na.rm = T),
            .groups = "drop")

results_to_plot <- results_to_plot %>%
  mutate(color_group = factor(
    case_when(
      bacteria_name == "S_aureus" & varying_param == 0.3 ~ 1,
      bacteria_name == "S_aureus" & varying_param == 0.6 ~ 2,
      bacteria_name == "S_aureus" & varying_param == 0.9 ~ 3,
      bacteria_name == "E_coli" & varying_param == 0.3 ~ 4,
      bacteria_name == "E_coli" & varying_param == 0.6 ~ 5,
      bacteria_name == "E_coli" & varying_param == 0.9 ~ 6
    ),
    labels = c(
      "0.3 - Bacteria 1",
      "0.6 - Bacteria 1",
      "0.9 - Bacteria 1",
      "0.3 - Bacteria 2",
      "0.6 - Bacteria 2",
      "0.9 - Bacteria 2"
    )
  ))

library(ggh4x)
ggplot(results_to_plot, aes(x = res_1y_wv_time, y = median, color = color_group))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = color_group), alpha = 0.2)+
  scale_color_manual(values = palette_combined_named)+
  scale_fill_manual(values = palette_combined_named)+
  facet_nested(bacteria_name+metric_name ~ name_renamed, scales = "free_y")+
  theme_bw()


ggplot(results_to_plot%>% filter(bacteria_name == "S_aureus", metric_name %in% c("resistance_carriage", "resistance_frequency")), aes(x = res_1y_wv_time, y = median, color = color_group))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = color_group), alpha = 0.2)+
  scale_color_manual(values = palette_combined_named)+
  scale_fill_manual(values = palette_combined_named)+
  facet_nested(bacteria_name+metric_name ~ name_renamed, scales = "free_y")+
  theme_bw()

ggplot(results_to_plot%>% filter(bacteria_name == "E_coli", metric_name %in% c("resistance_carriage", "resistance_frequency")), aes(x = res_1y_wv_time, y = median, color = color_group))+
  geom_line(linewidth = 1)+
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = color_group), alpha = 0.2)+
  scale_color_manual(values = palette_combined_named)+
  scale_fill_manual(values = palette_combined_named)+
  facet_nested(bacteria_name+metric_name ~ name_renamed, scales = "free_y")+
  theme_bw()
