library(lhs)
library(ggplot2)
library(patchwork)
library(dplyr)
library(deSolve)
library(purrr)
library(furrr)
library(tidyr)

source(here::here("scripts","utils","utils_lhs_and_prcc.R"))
source(here::here("scripts","utils","functions.R"))
source(here::here("scripts","model.R"))

n_cores = parallel::detectCores() - 1
plan(multisession, workers = n_cores)

number_of_samples <- 5000
set.seed(1234)

#### Define parameter distributions ####

# list of parameter distributions
param_distributions = list("betaC" = function(p) qunif(p, min = 0.01, max = 0.1),
                           "betaI" = function(p) qunif(p, min = 0.01, max = 0.1),
                           "f"     = function(p) qunif(p, min = 0.9, max = 1),
                           "as"    = function(p) qunif(p, min = 5*10^(-6), max = 10^(-4)),
                           "time_until_recovery_without_ATB_s" = function(p) qunif(p, min = 10, 50),
                           "dps"   = function(p) qunif(p, min = 30, max = 130),
                           "prob_bystander_exposure"             = function(p) qunif(p, min = 0.001, max = 0.01),
                           "time_until_decolo_by_bystander_ATB"  = function(p) qunif(p, min = 10, max = 30),
                           "prob_minority_strain_when_colonised" = function(p) qunif(p, min = 0, max = 0.30),
                           "prob_specific_exposure"              = function(p) qunif(p, min = 0.5, max = 1),
                           "time_until_decolo_by_specific_ATB"   = function(p) qunif(p, min = 3, max = 15),
                           "prob_minority_strain_when_infected" = function(p) qunif(p, min = 0, max = 0.30),
                           "prob_specific_exposure_r"              = function(p) qunif(p, min = 0.5, max = 1),
                           "time_until_decolo_by_specific_ATB_r"   = function(p) qunif(p, min = 5, max = 40))

params_fixed_values = data.frame(
  eps = 1,
  thetasr = 0,
  thetars = 0)

#### Create LHS dataframe ####


lhs_df <- generate_lhs(number_of_samples, param_distributions, params_fixed_values)

lhs_df <- set_resistant_and_sensitive_param_equal(lhs_df, c("dpr", "ar", "time_until_recovery_without_ATB_r"))

lhs_df_with_coexistence <- filter_coexistence_condition(lhs_df)

eq_results <- compute_equilibrium(lhs_df_with_coexistence)

# Simulate 1 year without vaccine for each parameter set at equilibrium
results_1y_wov <- eq_results %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov") %>%
  unnest_wider(res_1y_wov)


inputs <- results_1y_wov %>%
  select(all_of(names(param_distributions)))

outputs <- results_1y_wov %>%
  mutate(inccumIs = inccumIsnv,
         inccumIr = inccumIrnv,
         inccumI = inccumIs + inccumIr,
         prop_inccumIr = inccumIr/inccumI,
         prevCs = Csnv,
         prevCr = Crnv,
         prevC = prevCs + prevCr,
         prop_prevCr = prevCr/prevC)%>%
  select(c("inccumI", "inccumIs", "inccumIr","prop_inccumIr","prevC","prevCs","prevCr","prop_prevCr"))


#### Compute prcc with epiR ####

epi_prcc_df <- compute_prcc_epi_multi_outputs(inputs, outputs, c("inccumI", "inccumIs", "inccumIr","prop_inccumIr","prevC","prevCs","prevCr","prop_prevCr") )

save(epi_prcc_df, file = here::here("files",paste0("prcc_for_eq_",number_of_samples,".data")))

load(file = here::here("files",paste0("prcc_for_eq_",number_of_samples,".data")))

p_epi <- plot_prcc_epi_multi_color(epi_prcc_df, 4)
p_epi
ggsave(paste0("figures/figure2.png"), plot = p_epi, width = 14, height = 9)

# The code below allows you to have some parts of same graphic separately
# Only with the cumulative incidence of all infections
p_epi <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df, 4, 
                                                       cols_choice = c("inccumI"))+
  theme(legend.key.spacing.y = unit(0.5, "cm"))
p_epi
ggsave(paste0("figures/figure2_only_inccumI.png"), plot = p_epi, width = 7, height = 5)

# Only with the cumulative incidence of sensitive, resistant infections and proportion
p_epi <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df, 4, 
                                                       cols_choice = c("inccumIs","inccumIr", "prop_inccumIr"))+
  theme(legend.key.spacing.y = unit(0.5, "cm"))
p_epi
ggsave(paste0("figures/figure2_only_inccumIs_Ir_prop.png"), plot = p_epi, width = 14, height = 5)

# Run the code below if you want to have a graphic only with the cumulative incidence
p_epi <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df, 4, 
                                                       cols_choice = c("prop_inccumIr"))+
  theme(legend.key.spacing.y = unit(0.5, "cm"))
p_epi
ggsave(paste0("figures/figure2_only_prop_inccum.png"), plot = p_epi, width = 7, height = 5)

