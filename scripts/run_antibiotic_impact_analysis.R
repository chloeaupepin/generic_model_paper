library(purrr)
library(tidyr)
library(dplyr)
library(deSolve)
library(furrr)
library(stringr)

n_cores = parallel::detectCores() - 1
plan(multisession, workers = n_cores)


source(here::here("scripts", "model.R"))
source(here::here("scripts","utils", "functions.R"))

#### Define settings ####

## Population size
population_size = 100000

## Choose bacteria 
# S_aureus
# bacteria = "S_aureus"
# folder_name = "S_aureus"
# file_name  = "S_aureus_params10.csv"

# E_coli
bacteria = "E_coli"
folder_name = "E_coli"
file_name = "E_coli_params_primavera4.csv"

#### Load previous analysis equilibrium ####

load(here::here("files",folder_name,"equilibrium_results.RData"))

#### Simulate different antibiotic scenarios ####
# Simulate 1 year without vaccine for each parameter set at equilibrium

results_1y_wov <- eq_results %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov") %>%
  unnest_wider(res_1y_wov, names_sep = "_")

# Choose antibiotic scenario
antibiotic_scenarios_df <- tibble(prob_bystander_exposure_multiplier = c(0.9, 0.7, 0.5, 0.3,0.1)) %>%
  mutate(scenario_id = row_number())


# Add antibiotic parameters to each parameter set
df_for_antibiotic_simulations <- results_1y_wov %>%
  slice(rep(1:n(), each = nrow(antibiotic_scenarios_df))) %>%
  mutate(scenario_id = rep(1:nrow(antibiotic_scenarios_df), times = nrow(results_1y_wov))) %>%
  left_join(antibiotic_scenarios_df, by = "scenario_id")%>%
  mutate(prob_bystander_exposure_old = prob_bystander_exposure, .before = prob_bystander_exposure) %>%
  mutate(
    prob_bystander_exposure = prob_bystander_exposure * prob_bystander_exposure_multiplier
  ) %>%
  select(-prob_bystander_exposure_multiplier)

# Simulate one year with vaccine for each parameter set
results <- df_for_antibiotic_simulations %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov_antibiotic") %>%
  unnest_wider(res_1y_wov_antibiotic, names_sep = "_")


#### Results manipulation ####
# Compute resistance change
results_with_outputs_antibiotic <- results %>%
  mutate(res_1y_wov_inccumIrnv_sel = res_1y_wov_inccumIrnv + res_1y_wov_inccumSelectionOfResistantBySpecificAntibioForNonVaccinated + res_1y_wov_inccumSelectionOfResistantByBystanderAntibioForNonVaccinated,
         res_1y_wov_antibiotic_inccumIrnv_sel = res_1y_wov_antibiotic_inccumIrnv + res_1y_wov_antibiotic_inccumSelectionOfResistantBySpecificAntibioForNonVaccinated + res_1y_wov_antibiotic_inccumSelectionOfResistantByBystanderAntibioForNonVaccinated) %>%
  
  mutate(res_1y_wov_inccumI = res_1y_wov_inccumIrnv_sel+res_1y_wov_inccumIsnv,
         res_1y_wov_antibiotic_inccumI = res_1y_wov_antibiotic_inccumIrnv_sel+res_1y_wov_antibiotic_inccumIsnv) %>%
  
  mutate(res_1y_wov_prop_prevCr = (res_1y_wov_Crnv)/(res_1y_wov_Crnv + res_1y_wov_Csnv),
         res_1y_wov_antibiotic_prop_prevCr = (res_1y_wov_antibiotic_Crnv)/(res_1y_wov_antibiotic_Crnv + res_1y_wov_antibiotic_Csnv)) %>%
  
  mutate(res_1y_wov_prop_inccumIr = (res_1y_wov_inccumIrnv)/(res_1y_wov_inccumIrnv + res_1y_wov_inccumIsnv),
         res_1y_wov_antibiotic_prop_inccumIr = (res_1y_wov_antibiotic_inccumIrnv)/(res_1y_wov_antibiotic_inccumIrnv + res_1y_wov_antibiotic_inccumIsnv)) %>%
  
  mutate(prc_red_inccumI = prc_red(res_1y_wov_inccumI,res_1y_wov_antibiotic_inccumI),
         prc_red_inccumIr = prc_red(res_1y_wov_inccumIrnv_sel,res_1y_wov_antibiotic_inccumIrnv_sel),
         prc_red_inccumIs = prc_red(res_1y_wov_inccumIsnv,res_1y_wov_antibiotic_inccumIsnv),
         prc_red_prop_inccumIr = prc_red(res_1y_wov_prop_inccumIr, res_1y_wov_antibiotic_prop_inccumIr),
         prc_red_prevC = prc_red(res_1y_wov_Crnv + res_1y_wov_Csnv, res_1y_wov_antibiotic_Crnv + res_1y_wov_antibiotic_Csnv), 
         prc_red_prevCs = prc_red(res_1y_wov_Csnv, res_1y_wov_antibiotic_Csnv),
         prc_red_prevCr = prc_red(res_1y_wov_Crnv, res_1y_wov_antibiotic_Crnv ),
         prc_red_prop_prevCr = prc_red(res_1y_wov_prop_prevCr, res_1y_wov_antibiotic_prop_prevCr ))

save(results_with_outputs_antibiotic, file = here::here("files",folder_name,"results_with_outputs_antibiotic.RData"))

# Compute statistics
results_to_plot_antibiotic <- results_with_outputs_antibiotic %>%
  pivot_longer(cols = all_of(starts_with("prc_red")), names_to = "metric_name", values_to = "metric_value")%>%
  group_by(scenario_id, metric_name) %>%
  summarise(median = median(metric_value),
            q025 = quantile(metric_value, 0.025),
            q975 = quantile(metric_value, 0.975),
            .groups = "drop")


save(results_to_plot_antibiotic, file = here::here("files",folder_name,"results_to_plot_antibiotic.RData"))




