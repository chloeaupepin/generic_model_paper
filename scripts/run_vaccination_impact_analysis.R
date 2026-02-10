library(purrr)
library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2)
library(ggh4x)
library(furrr)
library(tidyr)
library(stringr)

n_cores = parallel::detectCores() - 1
plan(multisession, workers = n_cores)


source(here::here("scripts", "model.R"))
source(here::here("scripts","utils", "functions.R"))

set.seed(123)  # global seed for reproducibility

#### Define settings ####

## Choose bacteria 
# S_aureus
folder_name = "S_aureus"
file_name  = "S_aureus_params10.csv"

# E_coli
#folder_name = "E_coli"
#file_name = "E_coli_params_primavera2.csv"

## Choose number of bacteria to generate
n <- 50

## Is transmission by infected individuals allowed ? 
# For S_aureus no, for E_coli yes
# if transmission is allowed, it equals transmission by colonized individuals
transmission_by_infected = FALSE
#transmission_by_infected = TRUE


#### Define bacteria ####
Bacteria_params = read.csv(here::here("files",file_name)) 

# cols_to_noise include bacteria parameters and antibiotic parameters
# the relative fitness f isn't varied due to too much sensitivity of resistant proportion
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

df_bacteria <- add_variability(Bacteria_params, cols_to_noise, n ,0.05) 

df_bacteria <- set_resistant_and_sensitive_param_equal(df_bacteria, c("dpr", "ar", "time_until_recovery_without_ATB_r"))


if(transmission_by_infected){
  df_bacteria <- df_bacteria %>%
    mutate(betaI = betaC)
} else {
  df_bacteria <- df_bacteria %>%
    mutate(betaI = 0)
}


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
save(eq_results, file = here::here("files",folder_name,"equilibrium_results.RData"))

#### Simulation de un an sans vaccin et un an avec vaccin ####
# Simulate 1 year without vaccine for each parameter set at equilibrium

results_1y_wov <- eq_results %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov") %>%
  unnest_wider(res_1y_wov, names_sep = "_")

# Choose vaccination scenario
vaccine_scenarios_complete_df <- bind_rows(
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=0,"vfdr"=0,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0,name = "vftc"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vfis" = c(0.3,0.6,0.9)) %>% mutate("vfir" = vfis,"vftcs"=0, "vftcr"=0, "vfds"=0, "vfdr"=0, "vfrs"=0, "vfrr"=0, name = "vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfd"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfds"=0,"vfdr"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vfds" = c(0.3,0.6,0.9)) %>% mutate("vfdr" = vfds, "vfis"=vfds,"vfir"=vfds,"vftcs"=0, "vftcr"=0,"vfrs"=0, "vfrr"=0, name = "vfd_vfi"),
  expand.grid("Vperc" = c(0.1,0.3,0.5,0.7,0.9),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfrs"=0, "vfrr"=0, name = "vftc_vfd_vfi")
) %>% mutate(vaccine_id = row_number()) %>%
  mutate(vftis = 0, vftir = 0)

# Add vaccine parameters to each parameter set
df_for_vaccine_simulations <- add_vaccine_parameters(results_1y_wov, vaccine_scenarios_complete_df)

# Simulate one year with vaccine for each parameter set
results <- df_for_vaccine_simulations %>%
  apply_function_on_df(simulate_1y_with_vaccine, "res_1y_wv") %>%
  unnest_wider(res_1y_wv, names_sep = "_")


#### Results manipulation ####
# Compute resistance change
results_with_outputs <- compute_outputs(results)

results_with_outputs <- clean_vaccine_related_parameters(results_with_outputs)

# save results_with_outputs
save(results_with_outputs, file = here::here("files",folder_name,"results_with_outputs.RData"))

# load results_with_outputs
#load(here::here("files",folder_name,"results_with_outputs.RData"))

# Compute statistics
results_to_plot <- results_with_outputs %>%
  compute_statistics(starts_with(c("prc_red_", "averted_")))%>%
  mutate(
    population = case_when(
      str_detect(metric_name, "non_vaccinated") ~ "Non vaccinated",
      str_detect(metric_name, "_vaccinated") ~ "Vaccinated",
      TRUE ~ "Total population"
    )
  )

# save results_to_plot
save(results_to_plot, file = here::here("files",folder_name,"results_to_plot.RData"))
