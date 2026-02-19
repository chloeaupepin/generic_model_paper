library(lhs)
library(ggplot2)
library(patchwork)
library(dplyr)
library(deSolve)
library(purrr)
library(furrr)
library(tidyr)

n_cores = parallel::detectCores() - 1
plan(multisession, workers = n_cores)

source(here::here("scripts","utils", "utils_lhs_and_prcc.R"))
source(here::here("scripts","utils", "functions.R"))
source(here::here("scripts", "model.R"))

#### Define parameter distributions ####
param_distributions = list("betaC" = function(p) qunif(p, min = 0.01, max = 0.1),
                           "betaI" = function(p) qunif(p, min = 0.01, max = 0.1),
                           "f"     = function(p) qunif(p, min = 0.9, max = 1),
                           "as"    = function(p) qunif(p, min = 5*10^(-6), max = 10^(-4)),
                           "time_until_recovery_without_ATB_s" = function(p) qunif(p, min = 10, 50),
                           "dps"   = function(p) qunif(p, min = 30, max = 130),
                           "prob_bystander_exposure"             = function(p) qunif(p, min = 0.001, max = 0.01),
                           "time_until_decolo_by_bystander_ATB"  = function(p) qunif(p, min = 10, max = 30),
                           "prob_minority_strain_when_colonised" = function(p) qunif(p, min = 0, max = 0.5),
                           "prob_specific_exposure"              = function(p) qunif(p, min = 0.5, max = 1),
                           "time_until_decolo_by_specific_ATB"   = function(p) qunif(p, min = 3, max = 15),
                           "prob_minority_strain_when_infected" = function(p) qunif(p, min = 0, max = 0.5),
                           "prob_specific_exposure_r"              = function(p) qunif(p, min = 0.5, max = 1),
                           "time_until_decolo_by_specific_ATB_r"   = function(p) qunif(p, min = 5, max = 40),
                           "Vperc" = function(p) qunif(p, min = 0.1, max = 0.9),
                           "vftcs" = function(p) qunif(p, min = 0.3, max = 0.9),
                           "vfds" = function(p) qunif(p, min = 0.3, max = 0.9),
                           "vfis" = function(p) qunif(p, min = 0.3, max = 0.9),
                           "vfrs" = function(p) qunif(p, min = 0.3, max = 0.9)
)

params_fixed_values = data.frame(
  eps = 1,
  thetasr = 0,
  thetars = 0)


#### Create LHS dataframe ####
number_of_samples <- 5000
set.seed(1234)

lhs_df <- generate_lhs(number_of_samples, param_distributions, params_fixed_values)

lhs_df <- lhs_df %>% mutate(vftis = vftcs)

lhs_df <- set_resistant_and_sensitive_param_equal(lhs_df, c("dpr", "ar", "time_until_recovery_without_ATB_r",
                                                            "vftcr","vftir", "vfdr", "vfir", "vfrr"))

lhs_df_with_coexistence <- filter_coexistence_condition(lhs_df)

eq_results <- compute_equilibrium(lhs_df_with_coexistence)

# Simulate 1 year without vaccine 
results_1y_wov <- eq_results %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov") %>%
  unnest_wider(res_1y_wov, names_sep = "_")

# Simulate 1 year with vaccine 
results_1y_wov_wv <- results_1y_wov %>%
  apply_function_on_df(simulate_1y_with_vaccine, "res_1y_wv") %>%
  unnest_wider(res_1y_wv, names_sep = "_")


inputs <- results_1y_wov_wv %>%
  select(all_of(names(param_distributions)))

outputs <- results_1y_wov_wv %>%
  mutate(res_1y_wov_inccumIs = res_1y_wov_inccumIsnv,
         res_1y_wov_inccumIr = res_1y_wov_inccumIrnv,
         res_1y_wov_inccumI = res_1y_wov_inccumIs + res_1y_wov_inccumIr,
         res_1y_wov_prop_inccumIr = res_1y_wov_inccumIr/res_1y_wov_inccumI,
         res_1y_wov_prevCs = res_1y_wov_Csnv,
         res_1y_wov_prevCr = res_1y_wov_Crnv,
         res_1y_wov_prevC = res_1y_wov_prevCs + res_1y_wov_prevCr,
         res_1y_wov_prop_prevCr = res_1y_wov_prevCr/res_1y_wov_prevC,
         
         res_1y_wv_inccumIs = res_1y_wv_inccumIsnv +res_1y_wv_inccumIsv ,
         res_1y_wv_inccumIr = res_1y_wv_inccumIrnv+res_1y_wv_inccumIrv,
         res_1y_wv_inccumI = res_1y_wv_inccumIs + res_1y_wv_inccumIr,
         res_1y_wv_prop_inccumIr = res_1y_wv_inccumIr/res_1y_wv_inccumI,
         res_1y_wv_prevCs = res_1y_wv_Csnv+res_1y_wv_Csv,
         res_1y_wv_prevCr = res_1y_wv_Crnv + res_1y_wv_Crv,
         res_1y_wv_prevC = res_1y_wv_prevCs + res_1y_wv_prevCr,
         res_1y_wv_prop_prevCr = res_1y_wv_prevCr/res_1y_wv_prevC)%>%
  mutate(prc_red_inccumI = 100*(res_1y_wov_inccumI - res_1y_wv_inccumI)/res_1y_wov_inccumI,
         prc_red_inccumIs = 100*(res_1y_wov_inccumIs - res_1y_wv_inccumIs)/res_1y_wov_inccumIs,
         prc_red_inccumIr = 100*(res_1y_wov_inccumIr - res_1y_wv_inccumIr)/res_1y_wov_inccumIr,
         prc_red_prop_inccumIr = 100*(res_1y_wov_prop_inccumIr - res_1y_wv_prop_inccumIr)/res_1y_wov_prop_inccumIr,
         
         prc_red_prevC = 100*(res_1y_wov_prevC - res_1y_wv_prevC)/res_1y_wov_prevC,
         prc_red_prevCs = 100*(res_1y_wov_prevCs - res_1y_wv_prevCs)/res_1y_wov_prevCs,
         prc_red_prevCr = 100*(res_1y_wov_prevCr - res_1y_wv_prevCr)/res_1y_wov_prevCr,
         prc_red_prop_prevCr = 100*(res_1y_wov_prop_prevCr - res_1y_wv_prop_prevCr)/res_1y_wov_prop_prevCr
  ) %>%
  select(all_of(c("prc_red_inccumI", "prc_red_inccumIs", "prc_red_inccumIr","prc_red_prop_inccumIr",
                  "prc_red_prevC","prc_red_prevCs","prc_red_prevCr","prc_red_prop_prevCr")))

#### Compute prcc with epiR ####

epi_prcc_df <- compute_prcc_epi_multi_outputs(inputs, outputs, c("prc_red_inccumI", "prc_red_inccumIs", "prc_red_inccumIr","prc_red_prop_inccumIr",
                                                                 "prc_red_prevC","prc_red_prevCs","prc_red_prevCr","prc_red_prop_prevCr"
))

save(epi_prcc_df, file = here::here("files",paste0("prcc_for_vaccine_impact_",number_of_samples,".data")))

load(file = here::here("files",paste0("prcc_for_vaccine_impact_",number_of_samples,".data")))


p_epi <- plot_prcc_epi_multi_color(epi_prcc_df, 4)
p_epi
ggsave(paste0("figures/figure3.png"), plot = p_epi, width = 15, height = 10)

# Run the code below if you want to have a graphic only with the cumulative incidence
p_epi <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df, 4)+
  theme(legend.key.spacing.y = unit(0.5, "cm"))
p_epi
ggsave(paste0("figures/figure3_only_inccum.png"), plot = p_epi, width = 18, height = 6)

# Run the code below if you want to have a graphic only with the cumulative incidence
p_epi <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df, 4, cols_choice = c("prc_red_prevC","prc_red_prevCs","prc_red_prevCr","prc_red_prop_prevCr"))+
  theme(legend.key.spacing.y = unit(0.5, "cm"))
p_epi
ggsave(paste0("figures/figure3_only_prev.png"), plot = p_epi, width = 18, height = 6)
