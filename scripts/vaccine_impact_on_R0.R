
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(furrr)
library(ggrepel)

source(here::here("scripts","utils", "functions.R"))
source(here::here("scripts", "utils", "utils_R0.R"))
source(here::here("scripts", "utils", "plot_functions.R"))

n_cores = parallel::detectCores() - 1
plan(multisession, workers = n_cores)


## Population size
N = 100000

## Choose bacteria 
# Comment or uncomment line codes below and rerun everything for each bacteria 
# S_aureus
bacteria = "S_aureus"
folder_name = "S_aureus"
file_name  = "S_aureus_params_eu.csv"
transmission_by_infected = F

# E_coli
# bacteria = "E_coli"
# folder_name = "E_coli"
# file_name = "E_coli_params_primavera4.csv"
# transmission_by_infected = T

#### Load previous analysis equilibrium ####

load(here::here("files",folder_name,"equilibrium_results.RData"))

#### Add vaccination scenario ####
# Choose vaccination scenario
vaccine_scenarios_complete_df <- bind_rows(
  expand.grid("Vperc" = c(0,0.1,0.3,0.5,0.7,0.9,1),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=0,"vfdr"=0,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0,name = "vftc"),
  expand.grid("Vperc" = c(0,0.1,0.3,0.5,0.7,0.9,1),"vfds" = c(0.3,0.6,0.9)) %>% mutate("vfdr" = vfds, "vftcs"=0, "vftcr"=0,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0,name = "vfd"),
  expand.grid("Vperc" = c(0,0.1,0.3,0.5,0.7,0.9,1),"vfis" = c(0.3,0.6,0.9)) %>% mutate("vfir" = vfis,"vftcs"=0, "vftcr"=0, "vfds"=0, "vfdr"=0, "vfrs"=0, "vfrr"=0, name = "vfi"),
  expand.grid("Vperc" = c(0,0.1,0.3,0.5,0.7,0.9,1),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs,"vfis"=0,"vfir"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfd"),
  expand.grid("Vperc" = c(0,0.1,0.3,0.5,0.7,0.9,1),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfds"=0,"vfdr"=0,"vfrs" = 0, "vfrr" = 0, name = "vftc_vfi"),
  expand.grid("Vperc" = c(0,0.1,0.3,0.5,0.7,0.9,1),"vfds" = c(0.3,0.6,0.9)) %>% mutate("vfdr" = vfds, "vfis"=vfds,"vfir"=vfds,"vftcs"=0, "vftcr"=0,"vfrs"=0, "vfrr"=0, name = "vfd_vfi"),
  expand.grid("Vperc" = c(0,0.1,0.3,0.5,0.7,0.9,1),"vftcs" = c(0.3,0.6,0.9)) %>% mutate("vftcr" = vftcs, "vfds"=vftcs,"vfdr"=vftcs, "vfis"=vftcs,"vfir"=vftcs,"vfrs"=0, "vfrr"=0, name = "vftc_vfd_vfi")
) %>% mutate(vaccine_id = row_number()) 

if(transmission_by_infected){
  vaccine_scenarios_complete_df <- vaccine_scenarios_complete_df %>%
    mutate(vftis = vftcs, vftir = vftcr)
} else {
  vaccine_scenarios_complete_df <- vaccine_scenarios_complete_df %>%
    mutate(vftis = 0, vftir = 0)
}

df <- add_vaccine_parameters(eq_results, vaccine_scenarios_complete_df)

# Simulate one year with vaccine for each parameter set
results <- df %>%
  apply_function_on_df(SCISsrV.reprod_nbrs, "R0") %>%
  unnest_wider(R0) %>%
  mutate(R0 = pmax(R0s, R0r)) 

results <- results %>%
  mutate(R0s_R0r = R0s - R0r)


results <- clean_vaccine_related_parameters(results)

# Find best vaccine coverage
results_best_Vperc <- results %>%
  #filter(name == "vftc_vfd_vfi", vftcs == 0.9,  bacteria_id == 0, Vperc == 0) %>%
  apply_function_on_df(vaccine_coverage_threshold_for_R0, "best_Vperc")%>%
  mutate(best_Vperc = unlist(best_Vperc))

results_best_Vperc_stats <- results_best_Vperc %>%
  group_by(name_renamed, varying_param) %>%
  summarise(median = median(best_Vperc),
            q025 = quantile(best_Vperc, 0.025, na.rm = TRUE),
            q975 = quantile(best_Vperc, 0.975, na.rm = TRUE),
            .groups = "drop")

# Compute statistics
results_stats <- compute_statistics(results, c("R0", "R0s", "R0r", "R0s_R0r"))

#### Plots ####
chosen_palette = colors_and_shapes[[folder_name]]$palette
chosen_shape = colors_and_shapes[[folder_name]]$shape


plot_R0(data = results_stats,
        metric_name_to_plot = c("R0"),
        y_label = "Basic reproduction number",
        chosen_shape = chosen_shape,
        chosen_palette = chosen_palette,
        hrefs = c(1),
        vrefs = c()) +
  geom_point(data = results_best_Vperc_stats %>% filter(median > 0, median < 1), mapping = aes(x = median*100, y = 1, color = varying_param), shape = 8,size = 2, show.legend = FALSE)+
  geom_label_repel(data = results_best_Vperc_stats %>% filter(median > 0, median < 1), mapping = aes(x = median*100, y = 1, label = paste0(round(median*100),"%"), color = varying_param), size = 4, show.legend = F)
ggsave(here::here("figures",folder_name,"R0.png"), width = 14, height = 4)

plot_R0(data = results_stats,
        metric_name_to_plot = c("R0s_R0r"),
        y_label = "Difference between\nbasic reproduction number\nof sensitive and resistant",
        chosen_shape = chosen_shape,
        chosen_palette = chosen_palette,
        hrefs = c(0),
        vrefs = c())
ggsave(here::here("figures",folder_name,"diff_R0s_R0r.png"), width = 12, height = 5)


