library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)


source(here::here("scripts", "utils", "plot_functions.R"))

#### Load datasets ####
load(here::here("files","S_aureus","equilibrium_results.RData"))
eq_results_Saureus <- eq_results

load(here::here("files","E_coli","equilibrium_results.RData"))
eq_results_Ecoli <- eq_results

load(here::here("files", "S_aureus", "results_to_plot.RData"))
results_to_plot_Saureus <- results_to_plot

load(here::here("files", "E_coli", "results_to_plot.RData"))
results_to_plot_Ecoli <- results_to_plot

load(here::here("files", "S_aureus", "results_to_plot_antibiotic.RData"))
results_to_plot_antibiotic_Saureus <- results_to_plot_antibiotic

load(here::here("files", "E_coli", "results_to_plot_antibiotic.RData"))
results_to_plot_antibiotic_Ecoli <- results_to_plot_antibiotic

#### Plot equilibrium results description ####
# Saureus like
plot_histogram_eq_results(eq_results_Saureus,c(0.2,0.4),c(0,0.15) )
# Ecoli like
plot_histogram_eq_results(eq_results_Ecoli,c(0.9,1),c(0,0.2) )

plot_histogram_eq_results_two_bacteria(eq_results_Saureus, eq_results_Ecoli)


#### Plot vaccination impact analysis ####

folder_names = c("S_aureus", "E_coli")
list_of_results_to_plot = list(
  "S_aureus" = results_to_plot_Saureus, "E_coli" = results_to_plot_Ecoli
)


for(folder_name in folder_names){
  chosen_palette = colors_and_shapes[[folder_name]]$palette
  chosen_shape = colors_and_shapes[[folder_name]]$shape
  chosen_results_to_plot = list_of_results_to_plot[[folder_name]]
  

  plot_vaccine_metric(data = chosen_results_to_plot,
                      metric_name_to_plot = c("prc_red_inccumI", 
                                              "prc_red_inccumI_non_vaccinated", 
                                              "prc_red_inccumI_vaccinated"),
                      y_label = "Relative change due to vaccination (%)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0, -50),
                      vrefs = c(30),
                      ylim_values = c(-100, 1))
  ggsave(here::here("figures",folder_name,"prc_red_inccumI.png"), width = 12, height = 7)
  

  plot_vaccine_metric(data = chosen_results_to_plot,
                      metric_name_to_plot = c("prc_red_inccumIr", 
                                              "prc_red_inccumIr_non_vaccinated", 
                                              "prc_red_inccumIr_vaccinated"),
                      y_label = "Relative change due to vaccination (%)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0, -50),
                      vrefs = c(30),
                      ylim_values = c(-100, 1))
  ggsave(here::here("figures",folder_name,"prc_red_inccumIr.png"), width = 12, height = 7)
  
  plot_vaccine_metric(data = chosen_results_to_plot,
                      metric_name_to_plot = c("prc_red_inccumIs", 
                                              "prc_red_inccumIs_non_vaccinated", 
                                              "prc_red_inccumIs_vaccinated"),
                      y_label = "Relative change due to vaccination (%)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0, -50),
                      vrefs = c(30),
                      ylim_values = c(-100, 1))
  ggsave(here::here("figures",folder_name,"prc_red_inccumIs.png"), width = 12, height = 7)
  
  plot_vaccine_metric(data = chosen_results_to_plot,
                      metric_name_to_plot = c("prc_red_prevC", 
                                              "prc_red_prevC_non_vaccinated", 
                                              "prc_red_prevC_vaccinated"),
                      y_label = "Relative change due to vaccination (%)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0, -50),
                      vrefs = c(30))
  ggsave(here::here("figures",folder_name,"prc_red_prevC.png"), width = 12, height = 7)
  
  plot_vaccine_metric(data = chosen_results_to_plot,
                      metric_name_to_plot = c("prc_red_prevCr", 
                                              "prc_red_prevCr_non_vaccinated", 
                                              "prc_red_prevCr_vaccinated"),
                      y_label = "Relative change due to vaccination (%)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0, -50),
                      vrefs = c(30))
  ggsave(here::here("figures",folder_name,"prc_red_prevCr.png"), width = 12, height = 7)
  

  plot_vaccine_metric(data = chosen_results_to_plot,
                      metric_name_to_plot = c("prc_red_prop_prevCIr", 
                                              "prc_red_prop_prevCIr_non_vaccinated", 
                                              "prc_red_prop_prevCIr_vaccinated"),
                      y_label = "Relative change due to vaccination (%)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0),
                      vrefs = c(30))

  ggsave(here::here("figures",folder_name,"prc_red_prop_prevCIr.png"), width = 12, height = 7)
  

  plot_vaccine_metric(data = chosen_results_to_plot,
                      metric_name_to_plot = c("prc_red_prop_prevCr", 
                                              "prc_red_prop_prevCr_non_vaccinated", 
                                              "prc_red_prop_prevCr_vaccinated"),
                      y_label = "Relative change due to vaccination (%)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0),
                      vrefs = c(30))
  ggsave(here::here("figures",folder_name,paste0("prc_red_prop_prevCr.png")), width = 12, height = 7)
  

  plot_vaccine_metric(data = results_to_plot,
                      metric_name_to_plot = c("averted_inccumI"),
                      y_label = "Number of new infections averted in a year\n(per 100 000 individuals)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0),
                      vrefs = c(30))
  ggsave(here::here("figures",folder_name,"averted_inccumI.png"), width = 12, height = 5)
  
  plot_vaccine_metric(data = results_to_plot,
                      metric_name_to_plot = c("averted_inccumIr"),
                      y_label = "Number of new infections averted in a year\n(per 100 000 individuals)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0),
                      vrefs = c(30))
  ggsave(here::here("figures",folder_name,"averted_inccumIr.png"), width = 12, height = 5)
  
  plot_vaccine_metric(data = results_to_plot,
                      metric_name_to_plot = c("averted_inccumIs"),
                      y_label = "Number of new infections averted in a year\n(per 100 000 individuals)",
                      chosen_shape = chosen_shape,
                      chosen_palette = chosen_palette,
                      hrefs = c(0),
                      vrefs = c(30))
  ggsave(here::here("figures",folder_name,"averted_inccumIs.png"), width = 12, height = 5)
  
}

#### Plot combined vaccine impact ####

plot_vaccine_metric_both_bacteria(data1 = results_to_plot_Saureus,
                                  data2 = results_to_plot_Ecoli,
                                  metric_name_to_plot = c("averted_inccumI"),
                                  y_label = "Number of new infections averted in a year\n(per 100 000 individuals)",
                                  chosen_shapes = shapes,
                                  chosen_palette = palette_combined_named,
                                  hrefs = c(),
                                  vrefs = c(30),
                                  logscale = T)
ggsave(here::here("figures","figure5.png"), width = 14, height = 5)


#### Plot combined antibiotic impact ####

plot_antibiotic_metric_both_bacteria(data1 = results_to_plot_antibiotic_Saureus,
                                     data2 = results_to_plot_antibiotic_Ecoli,
                                     metric_name_to_plot = c("prc_red_inccumI", "prc_red_prop_prevCr"),
                                     chosen_shapes = shapes)
ggsave(here::here("figures","figure7.png"), width = 10, height = 5)


plot_antibiotic_metric_both_bacteria(data1 = results_to_plot_antibiotic_Saureus,
                                     data2 = results_to_plot_antibiotic_Ecoli,
                                     metric_name_to_plot = c("prc_red_inccumIr","prc_red_prevCr"),
                                     chosen_shapes = shapes)

plot_antibiotic_metric_both_bacteria(data1 = results_to_plot_antibiotic_Saureus,
                                     data2 = results_to_plot_antibiotic_Ecoli,
                                     metric_name_to_plot = c("prc_red_inccumI","prc_red_inccumIs","prc_red_inccumIr", "prc_red_prop_inccumIr",
                                                             "prc_red_prevC", "prc_red_prevCr", "prc_red_prevCs", "prc_red_prop_prevCr"),
                                     chosen_shapes = shapes,
                                     facet_cols = 4)

ggsave(paste0("figures/figure7.png"), width = 12, height = 7)


#### Create figure files ####

wd = getwd()

# Create figure 4 files
file.copy(from = file.path(wd, "figures", "S_aureus", "prc_red_inccumI.png"),
          to = file.path(wd, "figures", "figure4_Saureus.png"))

file.copy(from = file.path(wd, "figures", "E_coli", "prc_red_inccumI.png"),
          to = file.path(wd, "figures", "figure4_E_coli.png"))

# Create figure 6 files
file.copy(from = file.path(wd, "figures", "S_aureus", "prc_red_prop_prevCr.png"),
          to = file.path(wd, "figures", "figure6_Saureus.png"))

file.copy(from = file.path(wd, "figures", "E_coli", "prc_red_prop_prevCr.png"),
          to = file.path(wd, "figures", "figure6_E_coli.png"))

# Create figure S1 files
file.copy(from = file.path(wd, "figures", "S_aureus", "prc_red_inccumIr.png"),
          to = file.path(wd, "figures", "figureS1_Saureus.png"))

file.copy(from = file.path(wd, "figures", "E_coli", "prc_red_inccumIr.png"),
          to = file.path(wd, "figures", "figureS1_E_coli.png"))

