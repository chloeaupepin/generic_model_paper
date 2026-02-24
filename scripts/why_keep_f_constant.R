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


###############################################################################
# 1. S. aureus
###############################################################################


#### Define parameter distributions ####
Bacteria_params = read.csv(here::here("files","S_aureus_params10.csv")) 


# list of parameter distributions
dist <- function(p){
  qunif(p, min = 1-0.1, max = 1 + 0.1)
}
param_distributions = list("betaC" = function(p) Bacteria_params$betaC * dist(p),
                           "f"     = function(p) Bacteria_params$f * dist(p),
                           "as"    = function(p) Bacteria_params$as * dist(p),
                           "time_until_recovery_without_ATB_s" = function(p) Bacteria_params$time_until_recovery_without_ATB_s * dist(p),
                           "dps"   = function(p) Bacteria_params$dps * dist(p),
                           "prob_bystander_exposure"             = function(p) Bacteria_params$prob_bystander_exposure * dist(p),
                           "time_until_decolo_by_bystander_ATB"  = function(p) Bacteria_params$time_until_decolo_by_bystander_ATB * dist(p),
                           "prob_minority_strain_when_colonised" = function(p) Bacteria_params$prob_minority_strain_when_colonised * dist(p),
                           "prob_specific_exposure"              = function(p) Bacteria_params$prob_specific_exposure * dist(p),
                           "time_until_decolo_by_specific_ATB"   = function(p) Bacteria_params$time_until_decolo_by_specific_ATB * dist(p),
                           "prob_minority_strain_when_infected" = function(p) Bacteria_params$prob_minority_strain_when_infected * dist(p),
                           "prob_specific_exposure_r"              = function(p) Bacteria_params$prob_specific_exposure_r * dist(p),
                           "time_until_decolo_by_specific_ATB_r"   = function(p) Bacteria_params$time_until_decolo_by_specific_ATB_r * dist(p))

params_fixed_values = data.frame(
  betaI = 0,
  eps = 1,
  thetasr = 0,
  thetars = 0)

#### Create LHS dataframe ####

number_of_samples <- 5000
set.seed(1234)

lhs_df <- generate_lhs(number_of_samples, param_distributions, params_fixed_values)

lhs_df <- lhs_df %>%
  mutate(f = pmin(1, f ),
    prob_specific_exposure_r = pmin(1,prob_specific_exposure_r),
         prob_specific_exposure = pmin(1,prob_specific_exposure),
         prob_minority_strain_when_colonised = pmin(1,prob_minority_strain_when_colonised),
         prob_minority_strain_when_infected = pmin(1,prob_minority_strain_when_infected),
         prob_bystander_exposure = pmin(1,prob_bystander_exposure)
  )
  

lhs_df <- set_resistant_and_sensitive_param_equal(lhs_df, c("dpr", "ar", "time_until_recovery_without_ATB_r"))

#### Key parameters for the coexistence condition ####
lhs_df <- compute_antibiotic_associated_flux(lhs_df) %>%
  mutate(h = -f + (1/dpr + ar - etar * ar / (etar + gammaA_r)) /(1/dps + lambdaA + as - etas * as / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar/(etar + gammaA_r)))


inputs_condition <- lhs_df  %>%
  select(all_of(names(param_distributions)))

output_condition <- lhs_df %>% select(h)

epi_prcc_df_condition <- compute_prcc_epi(inputs_condition,output_condition) %>%
  rename(param = var, prcc = est) %>%
  mutate(output = "Coexistence condition", .before = 1)

#p_epi_condition <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df_condition , 1, c("Coexistence condition") )
#p_epi_condition


#### Key parameters for the equilibrium when condition is true ####

lhs_df_with_coexistence <- filter_coexistence_condition(lhs_df)

eq_results <- compute_equilibrium(lhs_df_with_coexistence)


# Simulate 1 year without vaccine for each parameter set at equilibrium
results_1y_wov <- eq_results %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov") %>%
  unnest_wider(res_1y_wov)


inputs <- results_1y_wov %>%
  select(all_of(names(param_distributions)))

outputs_of_interest = c("inccumI", "inccumIs", "inccumIr","prop_inccumIr",
                        "prevC","prevCs","prevCr","prop_prevCr")

outputs <- results_1y_wov %>%
  mutate(inccumIs = inccumIsnv,
         inccumIr = inccumIrnv,
         inccumI = inccumIs + inccumIr,
         prop_inccumIr = inccumIr/inccumI,
         prevCs = Csnv,
         prevCr = Crnv,
         prevC = prevCs + prevCr,
         prop_prevCr = prevCr/prevC)%>%
  select(all_of(outputs_of_interest))


# Compute prcc with epiR

epi_prcc_df <- compute_prcc_epi_multi_outputs(inputs, outputs, outputs_of_interest)

#p_epi <- plot_prcc_epi_multi_color(epi_prcc_df, 4)
#p_epi
#ggsave(paste0("figures/figureS1_2.png"), plot = p_epi, width = 14, height = 9)

#p_epi <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df , 1, c("prop_prevCr") )
#p_epi



#### Combined plot ####

epi_prcc_df_combined <- bind_rows(epi_prcc_df, epi_prcc_df_condition)

p_epi_combined <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df_combined , 2, c("Coexistence condition","prop_prevCr") ) +
  theme(legend.key.spacing.y = unit(0.5, "cm"))
p_epi_combined

ggsave(paste0("figures/why_keep_f_constant_Saureus.png"), plot = p_epi_combined, width = 11, height = 5)

###############################################################################
# 2. E. coli
###############################################################################


#### Define parameter distributions ####
Bacteria_params = read.csv(here::here("files","E_coli_params_primavera4.csv")) 

# List of columns to transform
cols_to_noise <- c(	
  "betaC","f","as","time_until_recovery_without_ATB_s","dps",
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

# Define distribution function
dist <- function(p){
  qunif(p, min = 1-0.1, max = 1 + 0.1)
}

# List of distributions
param_distributions <- map(cols_to_noise, ~ function(p) Bacteria_params[[.x]] * dist(p))
names(param_distributions) <- cols_to_noise

params_fixed_values = data.frame(
  eps = 1,
  thetasr = 0,
  thetars = 0)

#### Create LHS dataframe ####

number_of_samples <- 5000
set.seed(1234)

lhs_df <- generate_lhs(number_of_samples, param_distributions, params_fixed_values)

lhs_df <- lhs_df %>%
  mutate(f = pmin(1, f ),
         prob_specific_exposure_r = pmin(1,prob_specific_exposure_r),
         prob_specific_exposure = pmin(1,prob_specific_exposure),
         prob_minority_strain_when_colonised = pmin(1,prob_minority_strain_when_colonised),
         prob_minority_strain_when_infected = pmin(1,prob_minority_strain_when_infected),
         prob_bystander_exposure = pmin(1,prob_bystander_exposure)
  ) %>%
  mutate(betaI = betaC)


lhs_df <- set_resistant_and_sensitive_param_equal(lhs_df, c("dpr", "ar", "time_until_recovery_without_ATB_r"))

#### Key parameters for the coexistence condition ####
lhs_df <- compute_antibiotic_associated_flux(lhs_df, F) %>%
  mutate(h = -f + (1/dpr + ar - etar * ar / (etar + gammaA_r)) /(1/dps + lambdaA + as - etas * as / (etas + lambdaA_I + gammaA_s + psiA + phiA_I)+phiA)*(betaC + betaI*as/(etas + lambdaA_I + gammaA_s + psiA + phiA_I))/(betaC + betaI*ar/(etar + gammaA_r)))


inputs_condition <- lhs_df  %>%
  select(all_of(names(param_distributions)))

output_condition <- lhs_df %>% select(h)

epi_prcc_df_condition <- compute_prcc_epi(inputs_condition,output_condition) %>%
  rename(param = var, prcc = est) %>%
  mutate(output = "Coexistence condition", .before = 1)

#p_epi_condition <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df_condition , 1, c("Coexistence condition") )
#p_epi_condition


#### Key parameters for the equilibrium when condition is true ####

lhs_df_with_coexistence <- filter_coexistence_condition(lhs_df)

eq_results <- compute_equilibrium(lhs_df_with_coexistence)


# Simulate 1 year without vaccine for each parameter set at equilibrium
results_1y_wov <- eq_results %>%
  apply_function_on_df(simulate_1y_without_vaccine, "res_1y_wov") %>%
  unnest_wider(res_1y_wov)


inputs <- results_1y_wov %>%
  select(all_of(names(param_distributions)))

outputs_of_interest = c("inccumI", "inccumIs", "inccumIr","prop_inccumIr",
                        "prevC","prevCs","prevCr","prop_prevCr")

outputs <- results_1y_wov %>%
  mutate(inccumIs = inccumIsnv,
         inccumIr = inccumIrnv,
         inccumI = inccumIs + inccumIr,
         prop_inccumIr = inccumIr/inccumI,
         prevCs = Csnv,
         prevCr = Crnv,
         prevC = prevCs + prevCr,
         prop_prevCr = prevCr/prevC)%>%
  select(all_of(outputs_of_interest))


# Compute prcc with epiR

epi_prcc_df <- compute_prcc_epi_multi_outputs(inputs, outputs, outputs_of_interest)

#p_epi <- plot_prcc_epi_multi_color(epi_prcc_df, 4)
#p_epi
#ggsave(paste0("figures/figureS1_2.png"), plot = p_epi, width = 14, height = 9)

#p_epi <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df , 1, c("prop_prevCr") )
#p_epi



#### Combined plot ####

epi_prcc_df_combined <- bind_rows(epi_prcc_df, epi_prcc_df_condition)

p_epi_combined <- plot_prcc_epi_multi_color_only_inccum_or_prev(epi_prcc_df_combined , 2, c("Coexistence condition","prop_prevCr") )  +
  theme(legend.key.spacing.y = unit(0.5, "cm"))
p_epi_combined

ggsave(paste0("figures/why_keep_f_constant_Ecoli.png"), plot = p_epi_combined, width = 11, height = 5)

