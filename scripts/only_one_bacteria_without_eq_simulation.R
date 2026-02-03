library(purrr)
library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2)
library(ggh4x)

source(here::here("scripts", "model.R"))
source(here::here("scripts", "functions.R"))

#### Define bacteria ####
# S_aureus_params
file_name  = "S_aureus_params10.csv"
Bacteria_params = read.csv(here::here("files",file_name)) 

# E_coli_params
file_name = "E_coli_params_primavera1.csv"
Bacteria_params = read.csv(here::here("case_studies","data",file_name))

# Nombre de jeux de paramètres à générer
n <- 50

set.seed(123)  # global seed for reproducibility

cols_to_noise <- c(	
  "betaC",
  "as",
  "time_until_recovery_without_ATB_s",
  "dps",
  "prob_bystander_exposure",
  "time_until_decolo_by_bystander_ATB",
  "prob_minority_strain_when_colonised",
  "prob_specific_exposure",
  "time_until_decolo_by_specific_ATB",
  "prob_minority_strain_when_infected",
  "prob_specific_exposure_r",
  "time_until_decolo_by_specific_ATB_r"
)

df_bacteria <- add_variability(Bacteria_params, cols_to_noise, n ,0.05) 

df_bacteria <- set_resistant_and_sensitive_param_equal(df_bacteria, c("dpr", "ar", "time_until_recovery_without_ATB_r"))

df_bacteria <- df_bacteria %>%
  mutate(betaI = 0) #%>%
  #mutate(betaI = betaC)
  
df_bacteria <- df_bacteria %>%
  mutate(prob_specific_exposure_r = pmin(1,prob_specific_exposure_r),
         prob_specific_exposure = pmin(1,prob_specific_exposure),
         prob_minority_strain_when_colonised = pmin(1,prob_minority_strain_when_colonised),
         prob_minority_strain_when_infected = pmin(1,prob_minority_strain_when_infected),
         prob_bystander_exposure = pmin(1,prob_bystander_exposure)
  )

df_bacteria <- filter_coexistence_condition(df_bacteria)



#### Run analysis equilibrium ####

# Compute equilibrium with analytical expressions 
eq_results <- compute_equilibrium(df_bacteria)


#### Description des equilibres min max ####
View(eq_results %>%
       mutate(total_colonized = (eq_Cs + eq_Cr)/100000*100,
              resistance_ratio = eq_Cr / (eq_Cs + eq_Cr)*100,
              incidenceI = as*eq_Cs + ar*eq_Cr ))

# Compute max and min of total colonized and resistance ratio at equilibrium
View(eq_results %>%
       mutate(total_colonized = (eq_Cs + eq_Cr)/100000*100,
              resistance_ratio = eq_Cr / (eq_Cs + eq_Cr)*100,
              incidenceI = as*eq_Cs + ar*eq_Cr) %>%
       pivot_longer(cols = c(total_colonized, resistance_ratio, incidenceI), names_to = "metric", values_to = "value") %>%
       group_by(metric) %>%
       summarise(
         min = round(min(value),4),
         max = round(max(value),4),
       ))


#### Description des equilibres plots ####

# Saureus like
eq_results %>%
  mutate(total_colonized = (eq_Cs + eq_Cr)/100000,
         resistance_ratio = eq_Cr / (eq_Cs + eq_Cr)) %>%
  pivot_longer(cols = c(total_colonized, resistance_ratio), names_to = "metric", values_to = "value") %>%
  mutate(metric_name = case_when(
    metric == "total_colonized" ~ "Proportion colonized",
    metric == "resistance_ratio" ~ "Proportion of resistant\namong colonized"
  )) %>%
  mutate(metric_name = factor(metric_name, levels = c("Proportion colonized", "Proportion of resistant\namong colonized"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~metric_name, scales = "free_x")+
  labs(x = "Values",
       y = "Count") +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))+
  facetted_pos_scales(x = list(
    scale_x_continuous(limits = c(0.2,0.4)),
    scale_x_continuous(limits = c(0,0.15))
  ))
ggsave(here::here("exploration",file_name,"only_one_bacteria_equilibrium_descriptions.png"), width = 8, height = 4)

library(ggh4x)

# E coli like
eq_results %>%
  mutate(total_colonized = (eq_Cs + eq_Cr)/100000,
         resistance_ratio = eq_Cr / (eq_Cs + eq_Cr)) %>%
  pivot_longer(cols = c(total_colonized, resistance_ratio), names_to = "metric", values_to = "value") %>%
  mutate(metric_name = case_when(
    metric == "total_colonized" ~ "Proportion colonized",
    metric == "resistance_ratio" ~ "Proportion of resistant\namong colonized"
  )) %>%
  mutate(metric_name = factor(metric_name, levels = c("Proportion colonized", "Proportion of resistant\namong colonized"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~metric_name, scales = "free_x")+
  labs(x = "Values",
       y = "Count") +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))+
  facetted_pos_scales(x = list(
    scale_x_continuous(limits = c(0.9,1)),
    scale_x_continuous(limits = c(0,0.2))
  ))
ggsave(here::here("exploration",file_name,"only_one_bacteria_equilibrium_descriptions_Ecoli.png"), width = 8, height = 4)

# Both at the same time
eq_results_Saureus %>%
  mutate(bacteria_name = "Bacteria 1")%>%
  bind_rows(eq_results %>% mutate(bacteria_name = "Bacteria 2")) %>%
  mutate(total_colonized = (eq_Cs + eq_Cr)/100000,
         resistance_ratio = eq_Cr / (eq_Cs + eq_Cr)) %>%
  pivot_longer(cols = c(total_colonized, resistance_ratio), names_to = "metric", values_to = "value") %>%
  mutate(metric_name = case_when(
    metric == "total_colonized" ~ "Proportion colonized",
    metric == "resistance_ratio" ~ "Proportion of resistant\namong colonized"
  )) %>%
  mutate(metric_name = factor(metric_name, levels = c("Proportion colonized", "Proportion of resistant\namong colonized"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_grid(bacteria_name~metric_name, scales = "free_x")+
  labs(x = "Values",
       y = "Count") +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

#### Simulation de un an sans vaccin et un an avec vaccin ####
# Simulate 1 year without vaccine for each parameter set at equilibrium

eq_results <- eq_results %>%
  mutate(res_1y_wov = pmap(., function(...) {
    params <- list(...)      # On met tous les arguments de la ligne dans une liste
    simulate_1y_without_vaccine(params)      # On appelle la fonction avec cette liste
  })) %>%
  unnest_wider(res_1y_wov, names_sep = "_")

# Choose vaccination scenario
vaccine_scenarios_complete_df <- bind_rows(
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=0,"vfdr"=0,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0,name = "vftc"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vfis" = c(0.3,0.6,0.9)) %>% mutate("vfir" = vfis,"vftcs"=0, "vftcr"=0, "vfds"=0, "vfdr"=0, "vfrs"=0, "vfrr"=0, name = "vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfd"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfds"=0,"vfdr"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vfds" = c(0.3,0.6,0.9)) %>% mutate("vfdr" = vfds, "vfis"=vfds,"vfir"=vfds,"vftcs"=0, "vftcr"=0,"vfrs"=0, "vfrr"=0, name = "vfd_vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfrs"=0, "vfrr"=0, name = "vftc_vfd_vfi")
)


# Add vaccine parameters to each parameter set
df_for_vaccine <- eq_results %>%
  slice(rep(1:n(), each = nrow(vaccine_scenarios_complete_df))) %>%
  mutate(vaccine_id = rep(1:nrow(vaccine_scenarios_complete_df), times = nrow(eq_results))) %>%
  left_join(vaccine_scenarios_complete_df %>% mutate(vaccine_id = row_number()), by = "vaccine_id")

# Simulate one year with vaccine for each parameter set
results <- df_for_vaccine %>%
  mutate(res_1y_wv = pmap(., function(...) {
    params <- list(...)      # On met tous les arguments de la ligne dans une liste
    simulate_1y_with_vaccine(params)      # On appelle la fonction avec cette liste
  })) %>%
  unnest_wider(res_1y_wv, names_sep = "_")


#### MAnipule les résultats ####
# Compute resistance change
results <- results %>%
  mutate(res_1y_wov_inccumIrnv_tot = res_1y_wov_inccumIrnv + res_1y_wov_inccumSelectionOfResistantBySpecificAntibioForNonVaccinated + res_1y_wov_inccumSelectionOfResistantByBystanderAntibioForNonVaccinated,
         res_1y_wv_inccumIrnv_tot = res_1y_wv_inccumIrnv + res_1y_wv_inccumSelectionOfResistantBySpecificAntibioForNonVaccinated + res_1y_wv_inccumSelectionOfResistantByBystanderAntibioForNonVaccinated,
         res_1y_wv_inccumIrv_tot = res_1y_wv_inccumIrv + res_1y_wv_inccumSelectionOfResistantBySpecificAntibioForVaccinated + res_1y_wv_inccumSelectionOfResistantByBystanderAntibioForVaccinated) %>%
  mutate(resistance_change_prc = -(res_1y_wov_inccumIrnv_tot-(res_1y_wv_inccumIrnv_tot + res_1y_wv_inccumIrv_tot))/(res_1y_wov_inccumIrnv_tot)*100)%>%
  mutate(resistance_change_prc_non_vaccinated = -(res_1y_wov_inccumIrnv_tot*(1-Vperc)-res_1y_wv_inccumIrnv_tot)/(res_1y_wov_inccumIrnv_tot*(1-Vperc))*100)%>%
  mutate(resistance_change_prc_vaccinated = -(res_1y_wov_inccumIrnv_tot*Vperc-res_1y_wv_inccumIrv_tot)/(res_1y_wov_inccumIrnv_tot*Vperc)*100) %>%
  mutate(infections_change_prc = -(res_1y_wov_inccumIrnv_tot+res_1y_wov_inccumIsnv-(res_1y_wv_inccumIrnv_tot + res_1y_wv_inccumIrv_tot + res_1y_wv_inccumIsnv + res_1y_wv_inccumIsv))/(res_1y_wov_inccumIrnv_tot+res_1y_wov_inccumIsnv)*100)%>%
  mutate(infections_change_prc_non_vaccinated = -((res_1y_wov_inccumIrnv_tot+res_1y_wov_inccumIsnv)*(1-Vperc)-res_1y_wv_inccumIrnv_tot - res_1y_wv_inccumIsnv)/((res_1y_wov_inccumIrnv_tot+res_1y_wov_inccumIsnv)*(1-Vperc))*100)%>%
  mutate(infections_change_prc_vaccinated = -((res_1y_wov_inccumIrnv_tot+res_1y_wov_inccumIsnv)*Vperc-res_1y_wv_inccumIrv_tot-res_1y_wv_inccumIsv)/((res_1y_wov_inccumIrnv_tot+res_1y_wov_inccumIsnv)*Vperc)*100)

results<- results %>%
  mutate(ratio_rs_total_wov = (res_1y_wov_Crnv+res_1y_wov_Irnv)/(res_1y_wov_Crnv + res_1y_wov_Irnv + res_1y_wov_Csnv + res_1y_wov_Isnv),
         ratio_rs_total_wv = (res_1y_wv_Crnv+res_1y_wv_Irnv + res_1y_wv_Crv+res_1y_wv_Irv)/(res_1y_wv_Crnv + res_1y_wv_Irnv + res_1y_wv_Csnv + res_1y_wv_Isnv + res_1y_wv_Crv + res_1y_wv_Irv + res_1y_wv_Csv + res_1y_wv_Isv),
         ratio_rs_nv_wv = (res_1y_wv_Crnv+res_1y_wv_Irnv)/(res_1y_wv_Crnv + res_1y_wv_Irnv + res_1y_wv_Csnv + res_1y_wv_Isnv),
         ratio_rs_v_wv = (res_1y_wv_Crv+res_1y_wv_Irv)/(res_1y_wv_Crv + res_1y_wv_Irv + res_1y_wv_Csv + res_1y_wv_Isv)) %>%
  mutate(ratio_total_prc = 100*(ratio_rs_total_wv - ratio_rs_total_wov)/ratio_rs_total_wov,
         ratio_nv_prc = 100*(ratio_rs_nv_wv - ratio_rs_total_wov)/ratio_rs_total_wov,
         ratio_v_prc = 100*(ratio_rs_v_wv - ratio_rs_total_wov)/ratio_rs_total_wov)


results <- results %>%
  mutate(prev_resistance_change_prc = -(res_1y_wov_Crnv-(res_1y_wv_Crnv + res_1y_wv_Crv))/(res_1y_wov_Crnv)*100)%>%
  mutate(prev_resistance_change_prc_non_vaccinated = -(res_1y_wov_Crnv*(1-Vperc)-res_1y_wv_Crnv)/(res_1y_wov_Crnv*(1-Vperc))*100)%>%
  mutate(prev_resistance_change_prc_vaccinated = -(res_1y_wov_Crnv*Vperc-res_1y_wv_Crv)/(res_1y_wov_Crnv*Vperc)*100) %>%
  mutate(prev_total_change_prc = -(res_1y_wov_Crnv+res_1y_wov_Csnv-(res_1y_wv_Crnv + res_1y_wv_Crv + res_1y_wv_Csnv + res_1y_wv_Csv))/(res_1y_wov_Crnv+res_1y_wov_Csnv)*100)%>%
  mutate(prev_total_change_prc_non_vaccinated = -((res_1y_wov_Crnv+res_1y_wov_Csnv)*(1-Vperc)-res_1y_wv_Crnv - res_1y_wv_Csnv)/((res_1y_wov_Crnv+res_1y_wov_Csnv)*(1-Vperc))*100)%>%
  mutate(prev_total_change_prc_vaccinated = -((res_1y_wov_Crnv+res_1y_wov_Csnv)*Vperc-res_1y_wv_Crv-res_1y_wv_Csv)/((res_1y_wov_Crnv+res_1y_wov_Csnv)*Vperc)*100)


# Filter results
results_filtered <- results %>%
  #filter(R0s>1, R0r>1)%>%
  filter(equilibrium_outcome == "Coexistence") %>%
  mutate(varying_param = case_when(
    vftcs != 0 ~ vftcs,
    vfds != 0 ~ vfds,
    vfis != 0 ~ vfis,
    TRUE ~ 0
  )) %>% 
  filter(varying_param != 0) %>%
  mutate(name = factor(name, levels = c("vftc","vfi","vftc_vfd","vftc_vfi","vfd_vfi","vftc_vfd_vfi")))

results_filtered <- results_filtered %>%
  #mutate(R0 = pmax(R0s,R0r)) %>%
  mutate(Vperc_prc = Vperc*100) %>%
  mutate(Vperc_prc = factor(Vperc_prc, levels = c(10,30,50,70,90)))

# Rename vaccine names
results_filtered <- results_filtered %>%
  mutate(name_renamed= case_when(
    name == "vftc" ~ "va",
    name == "vfi" ~ "vi",
    name == "vftc_vfd" ~ "va_vd",
    name == "vftc_vfi" ~ "va_vi",
    name == "vfd_vfi" ~ "vd_vi",
    name == "vftc_vfd_vfi" ~ "va_vd_vi"
  )) %>%
  mutate(name_renamed = factor(name_renamed, levels = c("va","vi","va_vd","va_vi","vd_vi","va_vd_vi")))

# save results_to_plot
save(results_filtered, file = here::here("exploration",file_name,"only_one_bacteria_results_SAureus10.RData"))
save(results_filtered, file = here::here("exploration",file_name,"only_one_bacteria_results_EColi.RData"))

# load results_to_plot
#load(here::here("exploration",file_name,"only_one_bacteria_results_SAureus10.RData"))
#load(here::here("exploration",file_name,"only_one_bacteria_results_EColi.RData"))


#### Plot relative change results ####
results_to_plot <- results_filtered %>%
  pivot_longer(cols = c("resistance_change_prc", "resistance_change_prc_non_vaccinated","resistance_change_prc_vaccinated",
                        "infections_change_prc", "infections_change_prc_non_vaccinated","infections_change_prc_vaccinated",
                        "ratio_total_prc","ratio_nv_prc","ratio_v_prc",
                        "prev_resistance_change_prc", "prev_resistance_change_prc_non_vaccinated","prev_resistance_change_prc_vaccinated",
                        "prev_total_change_prc", "prev_total_change_prc_non_vaccinated","prev_total_change_prc_vaccinated"), 
               names_to = "metric",
               values_to = "metric_value")%>%
  mutate(metric_name = case_when(
    metric %in% c("resistance_change_prc", "resistance_change_prc_non_vaccinated","resistance_change_prc_vaccinated") ~ "Relative change of\nresistant infections",
    metric %in% c("infections_change_prc", "infections_change_prc_non_vaccinated","infections_change_prc_vaccinated") ~ "Relative change of\ntotal infections",
    metric %in% c("prev_resistance_change_prc", "prev_resistance_change_prc_non_vaccinated","prev_resistance_change_prc_vaccinated") ~ "Relative change of\nresistant colonisation",
    metric %in% c("prev_total_change_prc", "prev_total_change_prc_non_vaccinated","prev_total_change_prc_vaccinated") ~ "Relative change of\ntotal colonisation",
    metric %in% c("ratio_total_prc","ratio_nv_prc","ratio_v_prc") ~ "Relative change of\nresistant to sensitive prevalence ratio"
  ),
  population_target = case_when(
    metric %in% c("resistance_change_prc","infections_change_prc","prev_resistance_change_prc","prev_total_change_prc","ratio_total_prc") ~ "Total population",
    metric %in% c("resistance_change_prc_non_vaccinated","infections_change_prc_non_vaccinated","prev_resistance_change_prc_non_vaccinated","prev_total_change_prc_non_vaccinated","ratio_nv_prc") ~ "Non vaccinated",
    metric %in% c("resistance_change_prc_vaccinated","infections_change_prc_vaccinated","prev_resistance_change_prc_vaccinated","prev_total_change_prc_vaccinated","ratio_v_prc") ~ "Vaccinated",
  ))%>%
  mutate(metric_name = factor(metric_name, levels = c("Relative change of\ntotal infections","Relative change of\nresistant infections","Relative change of\ntotal colonisation","Relative change of\nresistant colonisation","Relative change of\nresistant to sensitive prevalence ratio")),
         population_target = factor(population_target, levels = c("Total population",  "Vaccinated","Non vaccinated")))%>%
  group_by(name_renamed, Vperc, varying_param, metric_name, population_target) %>%
  summarise(mean = mean(metric_value),
            min  = min(metric_value),
            max  = max(metric_value),
            sd = sd(metric_value),
            n = n(),
            .groups = "drop")%>%
  mutate(
   ci_low = mean - 1.96*sd,
   ci_high = mean + 1.96*sd
  )

# save results_to_plot
save(results_to_plot, file = here::here("exploration",file_name,"only_one_bacteria_results_to_plot_SAureus.RData"))
save(results_to_plot, file = here::here("exploration",file_name,"only_one_bacteria_results_to_plot_EColi.RData"))

palette_orange_dark <- c(
  "#FFB84D",  # moyen-clair
  "#FF9900",  # base
  "#CC7A00",  # moyen
  "#995C00",  # foncé
  "#663D00"   # très foncé
)

palette_orange_dark_3 <- palette_orange_dark[c(1,3,5)]

palette_pink_purple <- c(
  "#F3C6F3",  # très clair
  "#E19BE1",  # base
  "#C474C4",  # moyen
  "#9E4A9E",  # foncé
  "#733173"   # très foncé
)

palette_pink_purple_3 <- palette_pink_purple[c(1,3,5)]

# Choose palette for a given bacteria 
if(file_name == "S_aureus_params10.csv"){
  chosen_palette = palette_orange_dark_3
  chosen_shape = 21
} else{
  chosen_palette = palette_pink_purple_3
  chosen_shape = 23
}


plot_given_metric_population <- function(results_to_plot,metric_name_vector, population_target_vector,ymin = -100, ymax = 1){
  
  if(length(metric_name_vector)==1){
    p <- results_to_plot %>%
      filter(metric_name %in% metric_name_vector, population_target %in% population_target_vector) %>%
      ggplot(aes(x = Vperc*100)) +
      geom_ribbon(aes(ymin = min, ymax = max,
                      fill = factor(varying_param),
                      group = varying_param),
                  alpha = 0.2) +
      geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
      geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
      geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
      geom_line(aes(y = mean, color = factor(varying_param))) +
      geom_point(aes(y = mean, color = factor(varying_param)))+
      facet_nested(population_target ~ name_renamed)+
      scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
      scale_color_manual(values = palette_orange_dark_3)+
      scale_fill_manual(values = palette_orange_dark_3)+
      labs(x = "Vaccine coverage (%)",
           y = "Relative change due to vaccination (%)",
           color = "Vaccine efficacy",
           fill = "Vaccine efficacy")+
      theme_bw()+
      theme(legend.position = "bottom",
            strip.text = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.subtitle = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 14))
  }
  else if(length(population_target_vector)==1){
    p <- results_to_plot %>%
      filter(metric_name %in% metric_name_vector, population_target %in% population_target_vector) %>%
      ggplot(aes(x = Vperc*100)) +
      geom_ribbon(aes(ymin = min, ymax = max,
                      fill = factor(varying_param),
                      group = varying_param),
                  alpha = 0.2) +
      geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
      geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
      geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
      geom_line(aes(y = mean, color = factor(varying_param))) +
      geom_point(aes(y = mean, color = factor(varying_param)))+
      facet_nested(metric_name  ~ name_renamed)+
      scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
      scale_color_manual(values = palette_orange_dark_3)+
      scale_fill_manual(values = palette_orange_dark_3)+
      labs(x = "Vaccine coverage (%)",
           y = "Relative change due to vaccination (%)",
           color = "Vaccine efficacy",
           fill = "Vaccine efficacy")+
      theme_bw()+
      theme(legend.position = "bottom",
            strip.text = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.subtitle = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 14))
  }
  else{
    p <- results_to_plot %>%
      filter(metric_name %in% metric_name_vector, population_target %in% population_target_vector) %>%
      ggplot(aes(x = Vperc*100)) +
      geom_ribbon(aes(ymin = min, ymax = max,
                      fill = factor(varying_param),
                      group = varying_param),
                  alpha = 0.2) +
      geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
      geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
      geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
      geom_line(aes(y = mean, color = factor(varying_param))) +
      geom_point(aes(y = mean, color = factor(varying_param)))+
      facet_nested(metric_name + population_target ~ name_renamed)+
      scale_color_manual(values = palette_orange_dark_3)+
      scale_fill_manual(values = palette_orange_dark_3)+
      scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
      labs(x = "Vaccine coverage (%)",
           y = "Relative change due to vaccination (%)",
           color = "Vaccine efficacy",
           fill = "Vaccine efficacy")+
      theme_bw()+
      theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.subtitle = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 14))
  }
  
  return(p)
}

plot_given_metric_population(results_to_plot, c("Relative change of\ntotal infections","Relative change of\nresistant infections"), c("Total population", "Non vaccinated"))  

plot_given_metric_population(results_to_plot, c("Relative change of\ntotal infections"), c("Total population","Vaccinated"))  
ggsave(here::here("exploration","only_one_bacteria_infections_change_EColi.png"), width = 10, height = 5)

plot_given_metric_population(results_to_plot, c("Relative change of\nresistant to sensitive prevalence ratio"), c("Total population","Vaccinated"))  
ggsave(here::here("exploration","only_one_bacteria_ratio_change_EColi.png"), width = 10, height = 5)


plot_given_metric_population(results_to_plot, c("Relative change of\nresistant infections"), c("Total population","Non vaccinated"))


#### Exact plots for PACRI presentation 19/12 ####
p <- results_to_plot %>% filter(varying_param == 0.3) %>%
  filter(metric_name %in% c("Relative change of\ntotal infections"), 
         population_target %in% c("Total population","Vaccinated", "Non vaccinated")) %>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  #geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested(population_target ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Relative change due to vaccination (%)",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  ylim(-100,1)
p
ggsave(here::here("exploration",file_name,"only_one_bacteria_infections_change_one_efficacy.png"), plot = p, width = 10, height = 7)

p <- results_to_plot %>%
  filter(metric_name %in% c("Relative change of\ntotal infections"), 
         population_target %in% c("Total population","Vaccinated", "Non vaccinated")) %>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested(population_target ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Relative change due to vaccination (%)",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  ylim(-100,1)
p
ggsave(here::here("exploration",file_name,"only_one_bacteria_infections_change_all_efficacies.png"), width = 10, height = 7)

p <- results_to_plot %>%
  filter(metric_name %in% c("Relative change of\nresistant infections"), 
         population_target %in% c("Total population","Vaccinated", "Non vaccinated")) %>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested(population_target ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Relative change due to vaccination (%)",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  ylim(-100,1)
p
ggsave(here::here("exploration",file_name,"only_one_bacteria_resistance_change_all_efficacies.png"), width = 10, height = 7)


plot_given_metric_population(results_to_plot, c("Relative change of\nresistant to sensitive prevalence ratio"), c("Total population","Vaccinated","Non vaccinated"))  

p <- results_to_plot %>%
  filter(metric_name %in% c("Relative change of\nresistant to sensitive prevalence ratio"), 
         population_target %in% c("Total population","Vaccinated", "Non vaccinated")) %>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  #geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested(population_target ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Relative change due to vaccination (%)",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p

ggsave(here::here("exploration",file_name,"only_one_bacteria_ratio_change_all_pop.png"), width = 10, height = 7)


p <- results_filtered %>%
  mutate(inccum_infections_averted_total_s = res_1y_wov_inccumIsnv - res_1y_wv_inccumIsnv - res_1y_wv_inccumIsv,
         inccum_infections_averted_total_r = res_1y_wov_inccumIrnv_tot - res_1y_wv_inccumIrnv_tot - res_1y_wv_inccumIrv_tot,
         inccum_infections_averted_total_total = inccum_infections_averted_total_s + inccum_infections_averted_total_r) %>%
  group_by(name_renamed, Vperc, varying_param) %>%
  summarise(mean = mean(inccum_infections_averted_total_total),
            min = min(inccum_infections_averted_total_total),
            max = max(inccum_infections_averted_total_total),
            sd = sd(inccum_infections_averted_total_total),
            n = n()) %>%
  mutate(ci_low = mean - 1.96*sd,
         ci_high = mean + 1.96*sd)%>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  #geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested( ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Number of new infections averted in a year\n(per 100 000 individuals)",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p
ggsave(here::here("exploration",file_name,"only_one_bacteria_absolute_infections_change.png"), plot = p, width = 12, height = 5)

p <- results_to_plot %>%
  filter(metric_name %in% c("Relative change of\ntotal colonisation"), 
         population_target %in% c("Total population","Vaccinated", "Non vaccinated")) %>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested(population_target ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Relative change due to vaccination (%)",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  ylim(-100,1)
p
ggsave(here::here("exploration",file_name,"only_one_bacteria_prev_total_change_all_efficacies.png"), width = 10, height = 7)

p <- results_to_plot %>%
  filter(metric_name %in% c("Relative change of\nresistant colonisation"), 
         population_target %in% c("Total population","Vaccinated", "Non vaccinated")) %>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested(population_target ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Relative change due to vaccination (%)",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  ylim(-100,1)
p
ggsave(here::here("exploration",file_name,"only_one_bacteria_prev_resistant_change_all_efficacies.png"), width = 10, height = 7)

p <- results_filtered %>%
  mutate(prev_total = res_1y_wv_Crnv + res_1y_wv_Crv + res_1y_wv_Csnv + res_1y_wv_Csv) %>%
  group_by(name_renamed, Vperc, varying_param) %>%
  summarise(mean = mean(prev_total),
            min = min(prev_total),
            max = max(prev_total),
            sd = sd(prev_total),
            n = n()) %>%
  mutate(ci_low = mean - 1.96*sd,
         ci_high = mean + 1.96*sd)%>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  #geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested( ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Prevalence of total colonisation\nat the end of the year",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p
ggsave(here::here("exploration",file_name,"only_one_bacteria_absolute_prev_total.png"), plot = p, width = 10, height = 5)

p <- results_filtered %>%
  mutate(prev_resistance = res_1y_wv_Crnv + res_1y_wv_Crv) %>%
  group_by(name_renamed, Vperc, varying_param) %>%
  summarise(mean = mean(prev_resistance),
            min = min(prev_resistance),
            max = max(prev_resistance),
            sd = sd(prev_resistance),
            n = n()) %>%
  mutate(ci_low = mean-1.96*sd,
         ci_high = mean + 1.96*sd)%>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = factor(varying_param),
                  group = varying_param),
              alpha = 0.2) +
  #geom_hline(yintercept = -50, linetype="dashed", color = "grey")+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
  geom_point(aes(y = mean, color = factor(varying_param), fill = factor(varying_param)), shape = chosen_shape, size = 2)+
  geom_line(aes(y = mean, color = factor(varying_param))) +
  facet_nested( ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = chosen_palette)+
  scale_fill_manual(values = chosen_palette)+
  labs(x = "Vaccine coverage (%)",
       y = "Prevalence of resistant colonisation\nat the end of the year",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
p
ggsave(here::here("exploration",file_name,"only_one_bacteria_absolute_prev_resistance.png"), plot = p, width = 10, height = 5)



# Table of maximum reductions
get_table_max_reductions <- function(results_to_plot, Vperc_choice = 0.9, varying_param_choice = 0.9,metric_name_vector, population_target_vector){
  
  table_max <- results_to_plot %>%
    filter(metric_name %in% metric_name_vector, 
           population_target %in% population_target_vector,
           Vperc == Vperc_choice,
           varying_param == varying_param_choice) %>%
    group_by(name_renamed, metric_name, population_target) %>%
    summarise(max_reduction = min(mean),
              .groups = "drop")
  
  return(table_max)
}

get_table_max_reductions(results_to_plot, Vperc_choice = 0.9, varying_param_choice = 0.9,
                         metric_name_vector = c("Relative change of\nresistant infections"),
                         population_target_vector = c("Total population","Non vaccinated","Vaccinated"))


# Transform table to wide format with one column per name_renamed
table_max_wide <- get_table_max_reductions(results_to_plot, Vperc_choice = 0.9, varying_param_choice = 0.9,
                                            metric_name_vector = c("Relative change of\nresistant infections","Relative change of\ntotal infections","Relative change of\nresistant to sensitive prevalence ratio"),
                                            population_target_vector = c("Total population","Non vaccinated","Vaccinated")) %>%
  mutate(max_reduction = round(max_reduction, 1))%>%
  pivot_wider(names_from = name_renamed, values_from = max_reduction)
  

# Heatmap of max reductions
get_table_max_reductions(results_to_plot, Vperc_choice = 0.9, varying_param_choice = 0.9,
                         metric_name_vector = c("Relative change of\nresistant infections","Relative change of\ntotal infections","Relative change of\nresistant to sensitive prevalence ratio"),
                         population_target_vector = c("Total population","Non vaccinated","Vaccinated")) %>%
  ggplot(aes(x = name_renamed, y = population_target)) +
  geom_tile(aes(fill = max_reduction), color = "white") +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred", midpoint = 0,
                       name = "Max reduction (%)") +
  geom_text(aes(label = round(max_reduction, 1)), color = "black", size = 5) +
  labs(x = "Vaccine name", y = "") +
  facet_wrap(~ metric_name) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))


#### Combining both bacteria ####
file_name_1  = "S_aureus_params10.csv"
load(here::here("exploration",file_name_1,"only_one_bacteria_results_to_plot_SAureus.RData"))
results_to_plot_SAureus <- results_to_plot
file_name_2 = "E_coli_params_primavera1.csv"
load(here::here("exploration",file_name_2, "only_one_bacteria_results_to_plot_EColi.RData"))
results_to_plot_EColi <- results_to_plot

results_to_plot_both <- results_to_plot_SAureus %>%
  mutate(bacteria = "1") %>%
  bind_rows(results_to_plot_EColi %>% mutate(bacteria = "2"))

library(ggnewscale)



palette_combined_named = c(
  "0.3 - Bacteria 1" = palette_orange_dark[1],
  "0.6 - Bacteria 1" = palette_orange_dark[3],
  "0.9 - Bacteria 1" = palette_orange_dark[5],
  "0.3 - Bacteria 2" = palette_pink_purple[1],
  "0.6 - Bacteria 2" = palette_pink_purple[3],
  "0.9 - Bacteria 2" = palette_pink_purple[5]
)

plot_given_metric_population_both_bacteria <- function(results_to_plot_both, 
                                                       metric_name_filter, 
                                                       vaccine_efficacy_filter, 
                                                       ylim_min, ylim_max, 
                                                       hlines = c(-50, 0),
                                                       logscale=F){
  results_to_plot_both %>%
    filter(metric_name %in% metric_name_filter, 
           population_target %in% c("Total population","Vaccinated", "Non vaccinated")) %>%
    mutate(color_group = factor(
      case_when(
        bacteria == "1" & varying_param == 0.3 ~ 1,
        bacteria == "1" & varying_param == 0.6 ~ 2,
        bacteria == "1" & varying_param == 0.9 ~ 3,
        bacteria == "2" & varying_param == 0.3 ~ 4,
        bacteria == "2" & varying_param == 0.6 ~ 5,
        bacteria == "2" & varying_param == 0.9 ~ 6
      ),
      labels = c(
        "0.3 - Bacteria 1",
        "0.6 - Bacteria 1",
        "0.9 - Bacteria 1",
        "0.3 - Bacteria 2",
        "0.6 - Bacteria 2",
        "0.9 - Bacteria 2"
      )
    )) %>%
    mutate(color = case_when(
      bacteria == "1" & varying_param == 0.3 ~ palette_orange_dark_3[1],
      bacteria == "1" & varying_param == 0.6 ~ palette_orange_dark_3[3],
      bacteria == "1" & varying_param == 0.9 ~ palette_orange_dark_3[5],
      bacteria == "2" & varying_param == 0.3 ~ palette_pink_purple_3[1],
      bacteria == "2" & varying_param == 0.6 ~ palette_pink_purple_3[3],
      bacteria == "2" & varying_param == 0.9 ~ palette_pink_purple_3[5]
    )) %>%
    filter(varying_param %in% vaccine_efficacy_filter) %>%
    ggplot(aes(x = Vperc*100)) +
    geom_ribbon(aes(ymin = min, ymax = max,
                    fill = color_group,
                    group = color_group),
                alpha = 0.2) +
    geom_hline(yintercept = hlines[1], linetype="dashed", color = "grey")+
    geom_hline(yintercept = hlines[2], linetype="dashed", color = "grey")+
    geom_vline(xintercept = 30, linetype="dashed", color = "grey")+
    geom_line(aes(y = mean, color = color_group)) +
    geom_point(aes(y = mean, color = color_group))+
    facet_nested(population_target ~ name_renamed)+
    {if(logscale) scale_y_log10()}+
    scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
    scale_color_manual(values = palette_combined_named)+
    scale_fill_manual(values = palette_combined_named)+
    labs(x = "Vaccine coverage (%)",
         y = "Relative change due to vaccination (%)",
         color = "Vaccine efficacy",
         fill = "Vaccine efficacy")+
    theme_bw()+
    theme(legend.position = "bottom",
          strip.text = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.subtitle = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))+
    {if (ylim_min != ylim_max) ylim(ylim_min,ylim_max)}+
    guides(fill = guide_legend(nrow = 2, byrow = TRUE), color = guide_legend(nrow = 2, byrow = TRUE))
}

plot_given_metric_population_both_bacteria(results_to_plot_both,c("Relative change of\ntotal infections"),c(0.3,0.6,0.9), -100, 1 )
ggsave(here::here("exploration",file_name,"only_one_bacteria_infections_change_both_bacteria_all_coverage_all.png"), width = 10, height = 7)

plot_given_metric_population_both_bacteria(results_to_plot_both,c("Relative change of\ntotal infections"),c(0.3,0.6), -100, 1 )
ggsave(here::here("exploration",file_name,"only_one_bacteria_infections_change_both_bacteria_all_coverage_0.3_0.6.png"), width = 10, height = 7)

plot_given_metric_population_both_bacteria(results_to_plot_both,c("Relative change of\ntotal infections"),c(0.3), -100, 1 )
ggsave(here::here("exploration",file_name,"only_one_bacteria_infections_change_both_bacteria_all_coverage_0.3.png"), width = 10, height = 7)

plot_given_metric_population_both_bacteria(results_to_plot_both,
                                           c("Relative change of\nresistant to sensitive prevalence ratio"),
                                           c(0.3,0.6,0.9), 0, 0, hlines = c(-20, 0))
ggsave(here::here("exploration",file_name,"only_one_bacteria_ratio_change_both_bacteria_all_coverage_all.png"), width = 10, height = 7)


#### Combining both bacteria ####
file_name_1  = "S_aureus_params10.csv"
load(here::here("exploration",file_name_1,"only_one_bacteria_results_SAureus.RData"))
results_to_plot_SAureus <- results_filtered
file_name_2 = "E_coli_params_primavera1.csv"
load(here::here("exploration",file_name_2, "only_one_bacteria_results_EColi.RData"))
results_to_plot_EColi <- results_filtered

file_name = paste0(file_name_1,"_",file_name_2)

results_to_plot_both <- results_to_plot_SAureus %>%
  mutate(bacteria = "1") %>%
  bind_rows(results_to_plot_EColi %>% mutate(bacteria = "2")) %>%
  mutate(inccum_infections_averted_total_s = res_1y_wov_inccumIsnv - res_1y_wv_inccumIsnv - res_1y_wv_inccumIsv,
         inccum_infections_averted_total_r = res_1y_wov_inccumIrnv_tot - res_1y_wv_inccumIrnv_tot - res_1y_wv_inccumIrv_tot,
         inccum_infections_averted_total_total = inccum_infections_averted_total_s + inccum_infections_averted_total_r) 


results_to_plot_both %>%
  pivot_longer(cols = c("res_1y_wov_inccumIsnv","res_1y_wov_inccumIrnv"),
               names_to = "metric_name",
               values_to = "value")%>%
  group_by(bacteria, metric_name) %>%
  summarise(mean = mean(value),
            min = min(value),
            max = max(value),
            sd = sd(value),
            n = n()
  )%>%
  mutate(ci_low = mean - 1.96*sd,
         ci_high = mean + 1.96*sd)

p2 <- results_to_plot_both %>%
  pivot_longer(cols = c("inccum_infections_averted_total_s","inccum_infections_averted_total_r"),
               names_to = "metric_name",
               values_to = "metric_value"
  ) %>%
  mutate(metric_name = case_when(
    metric_name == "inccum_infections_averted_total_s" ~ "Sensitive infections",
    metric_name == "inccum_infections_averted_total_r" ~ "Resistant infections",
    
  ))%>%
  group_by(name_renamed, Vperc, varying_param, bacteria, metric_name) %>%
  summarise(mean = mean(metric_value),
            min = min(metric_value),
            max = max(metric_value)) %>%
  ungroup()%>%
  mutate(color_group = factor(
    case_when(
      bacteria == "1" & varying_param == 0.3 ~ 1,
      bacteria == "1" & varying_param == 0.6 ~ 2,
      bacteria == "1" & varying_param == 0.9 ~ 3,
      bacteria == "2" & varying_param == 0.3 ~ 4,
      bacteria == "2" & varying_param == 0.6 ~ 5,
      bacteria == "2" & varying_param == 0.9 ~ 6
    ),
    labels = c(
      "0.3 - Bacteria 1",
      "0.6 - Bacteria 1",
      "0.9 - Bacteria 1",
      "0.3 - Bacteria 2",
      "0.6 - Bacteria 2",
      "0.9 - Bacteria 2"
    )
  )) %>%
  mutate(color = case_when(
    bacteria == "1" & varying_param == 0.3 ~ palette_orange_dark_3[1],
    bacteria == "1" & varying_param == 0.6 ~ palette_orange_dark_3[3],
    bacteria == "1" & varying_param == 0.9 ~ palette_orange_dark_3[5],
    bacteria == "2" & varying_param == 0.3 ~ palette_pink_purple_3[1],
    bacteria == "2" & varying_param == 0.6 ~ palette_pink_purple_3[3],
    bacteria == "2" & varying_param == 0.9 ~ palette_pink_purple_3[5]
  )) %>%
  #filter(varying_param == 0.3) %>%
  #filter(varying_param %in% c(0.3,0.6)) %>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = min, ymax = max,
                  fill = color_group,
                  group = color_group),
              alpha = 0.2) +
  geom_line(aes(y = mean, color = color_group)) +
  geom_point(aes(y = mean, color = color_group))+
  facet_nested(metric_name ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = palette_combined_named)+
  scale_fill_manual(values = palette_combined_named)+
  scale_y_log10()+
  labs(x = "Vaccine coverage (%)",
       y = "Number of new infections averted in a year",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE), color = guide_legend(nrow = 2, byrow = TRUE))
p2

p2 <- results_to_plot_both %>%
  group_by(name_renamed, Vperc, varying_param, bacteria) %>%
  summarise(mean = mean(inccum_infections_averted_total_total),
            min = min(inccum_infections_averted_total_total),
            max = max(inccum_infections_averted_total_total),
            sd = sd(inccum_infections_averted_total_total)) %>%
  mutate(ci_low = mean - 1.96*sd,
         ci_high = mean + 1.96*sd)%>%
  ungroup()%>%
  mutate(color_group = factor(
    case_when(
      bacteria == "1" & varying_param == 0.3 ~ 1,
      bacteria == "1" & varying_param == 0.6 ~ 2,
      bacteria == "1" & varying_param == 0.9 ~ 3,
      bacteria == "2" & varying_param == 0.3 ~ 4,
      bacteria == "2" & varying_param == 0.6 ~ 5,
      bacteria == "2" & varying_param == 0.9 ~ 6
    ),
    labels = c(
      "0.3 - Bacteria 1",
      "0.6 - Bacteria 1",
      "0.9 - Bacteria 1",
      "0.3 - Bacteria 2",
      "0.6 - Bacteria 2",
      "0.9 - Bacteria 2"
    )
  )) %>%
  mutate(color = case_when(
    bacteria == "1" & varying_param == 0.3 ~ palette_orange_dark_3[1],
    bacteria == "1" & varying_param == 0.6 ~ palette_orange_dark_3[3],
    bacteria == "1" & varying_param == 0.9 ~ palette_orange_dark_3[5],
    bacteria == "2" & varying_param == 0.3 ~ palette_pink_purple_3[1],
    bacteria == "2" & varying_param == 0.6 ~ palette_pink_purple_3[3],
    bacteria == "2" & varying_param == 0.9 ~ palette_pink_purple_3[5]
  )) %>%
  ggplot(aes(x = Vperc*100)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high,
                  fill = color_group,
                  group = color_group),
              alpha = 0.2) +
  geom_point(aes(y = mean, color = color_group,fill = color_group, shape = color_group), size = 2)+
  geom_line(aes(y = mean, color = color_group)) +
  facet_nested( ~ name_renamed)+
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
  scale_color_manual(values = palette_combined_named)+
  scale_fill_manual(values = palette_combined_named)+
  scale_shape_manual(values = c(21,21,21,23,23, 23))+
  scale_y_log10()+
  labs(x = "Vaccine coverage (%)",
       y = "Number of new infections averted in a year\n(per 100 000 individuals)",
       color = "Vaccine efficacy",
       fill = "Vaccine efficacy")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  guides(fill = guide_legend(nrow = 2, 
                             byrow = TRUE,
                             override.aes = list(shape = ifelse(grepl("Bacteria 1", levels(p2$data$color_group)), 21, 23))), 
         color = guide_legend(nrow = 2, 
                              byrow = TRUE), 
         shape = "none")
p2
ggsave(here::here("exploration",file_name,"only_one_bacteria_absolute_infections_change_both_bacteria_all_coverage_all.png"), plot = p2, width = 10, height = 5)

