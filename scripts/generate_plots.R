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
# # Saureus like
# plot_histogram_eq_results(eq_results_Saureus,c(0.2,0.4),c(0,0.15) )
# 
# View(eq_results_Saureus %>%
#        mutate(total_colonized = (eq_Csnv + eq_Crnv)/100000*100,
#               resistance_ratio = eq_Crnv / (eq_Csnv + eq_Crnv)*100,
#               incidenceI = as*eq_Csnv + ar*eq_Crnv) %>%
#        pivot_longer(cols = c(total_colonized, resistance_ratio, incidenceI), names_to = "metric", values_to = "value") %>%
#        group_by(metric) %>%
#        summarise(
#          min = round(min(value),4),
#          max = round(max(value),4),
#        ))
# 
# # Ecoli like
# plot_histogram_eq_results(eq_results_Ecoli,c(0.9,1),c(0,0.2) )
# 
# View(eq_results_Ecoli %>%
#        mutate(total_colonized = (eq_Csnv + eq_Crnv)/100000*100,
#               resistance_ratio = eq_Crnv / (eq_Csnv + eq_Crnv)*100,
#               incidenceI = as*eq_Csnv + ar*eq_Crnv) %>%
#        pivot_longer(cols = c(total_colonized, resistance_ratio, incidenceI), names_to = "metric", values_to = "value") %>%
#        group_by(metric) %>%
#        summarise(
#          min = round(min(value),4),
#          max = round(max(value),4),
#        ))

# Combined
# plot_histogram_eq_results_two_bacteria(eq_results_Saureus, eq_results_Ecoli)




#### Plot vaccination impact analysis for S_aureus ####

folder_name = "S_aureus"
chosen_palette = colors_and_shapes[[folder_name]]$palette
chosen_shape = colors_and_shapes[[folder_name]]$shape
chosen_results_to_plot = results_to_plot_Saureus


plot_vaccine_metric(data = chosen_results_to_plot,
                    metric_name_to_plot = c("prc_red_inccumI", 
                                            "prc_red_inccumI_non_vaccinated", 
                                            "prc_red_inccumI_vaccinated"),
                    vaccine_effects_to_plot = c(0.3, 0.6, 0.9),
                    y_label = "Relative change due to vaccination (%)",
                    chosen_shape = chosen_shape,
                    chosen_palette = chosen_palette,
                    hrefs = c(0, -50),
                    vrefs = c(30),
                    ylim_values = c(-100, 1))
ggsave(here::here("figures",folder_name,"prc_red_inccumI.png"), width = 12, height = 7)
ggsave(here::here("figures","figure4_S_aureus.png"), width = 12, height = 7)

plot_vaccine_metric(data = chosen_results_to_plot,
                    metric_name_to_plot = c("prc_red_inccumI", 
                                            "prc_red_inccumI_non_vaccinated", 
                                            "prc_red_inccumI_vaccinated"),
                    vaccine_effects_to_plot = c(0.3),
                    y_label = "Relative change due to vaccination (%)",
                    chosen_shape = chosen_shape,
                    chosen_palette = chosen_palette,
                    hrefs = c(0, -50),
                    vrefs = c(30),
                    ylim_values = c(-100, 1))
ggsave(here::here("figures",folder_name,"prc_red_inccumI_1.png"), width = 12, height = 7)


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
ggsave(here::here("figures","figureS4_S_aureus.png"), width = 12, height = 7)

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
ggsave(here::here("figures","figureS5_S_aureus.png"), width = 12, height = 7)

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
ggsave(here::here("figures","figure6_S_aureus.png"), width = 12, height = 7)

plot_vaccine_metric(data = chosen_results_to_plot,
                    metric_name_to_plot = c("prc_red_prop_prevCr"),
                    y_label = "Relative change due to vaccination (%)",
                    chosen_shape = chosen_shape,
                    chosen_palette = chosen_palette,
                    hrefs = c(0),
                    vrefs = c(30))
ggsave(here::here("figures",folder_name,paste0("prc_red_prop_prevCr_total.png")), width = 14, height = 5)


plot_vaccine_metric(data = chosen_results_to_plot,
                    metric_name_to_plot = c("prc_red_prop_inccumIr", 
                                            "prc_red_prop_inccumIr_non_vaccinated", 
                                            "prc_red_prop_inccumIr_vaccinated"),
                    y_label = "Relative change due to vaccination (%)",
                    chosen_shape = chosen_shape,
                    chosen_palette = chosen_palette,
                    hrefs = c(0),
                    vrefs = c(30))
plot_vaccine_metric(data = chosen_results_to_plot,
                    metric_name_to_plot = c("prc_red_prop_inccumIr_sel", 
                                            "prc_red_prop_inccumIr_sel_non_vaccinated", 
                                            "prc_red_prop_inccumIr_sel_vaccinated"),
                    y_label = "Relative change due to vaccination (%)",
                    chosen_shape = chosen_shape,
                    chosen_palette = chosen_palette,
                    hrefs = c(0),
                    vrefs = c(30))
ggsave(here::here("figures",folder_name,paste0("prc_red_prop_inccumIr.png")), width = 12, height = 7)

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
  

#### Plot vaccination impact analysis for E_coli ####

folder_name = "E_coli"
chosen_palette = colors_and_shapes[[folder_name]]$palette
chosen_shape = colors_and_shapes[[folder_name]]$shape
chosen_results_to_plot = results_to_plot_Ecoli


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
ggsave(here::here("figures","figure4_E_coli.png"), width = 12, height = 7)

plot_vaccine_metric(data = chosen_results_to_plot,
                    metric_name_to_plot = c("prc_red_inccumI", 
                                            "prc_red_inccumI_non_vaccinated", 
                                            "prc_red_inccumI_vaccinated"),
                    vaccine_effects_to_plot = c(0.3),
                    y_label = "Relative change due to vaccination (%)",
                    chosen_shape = chosen_shape,
                    chosen_palette = chosen_palette,
                    hrefs = c(0, -50),
                    vrefs = c(30),
                    ylim_values = c(-100, 1))
ggsave(here::here("figures",folder_name,"prc_red_inccumI_1.png"), width = 12, height = 7)


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
ggsave(here::here("figures","figureS4_E_coli.png"), width = 12, height = 7)

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
ggsave(here::here("figures","figureS5_E_coli.png"), width = 12, height = 7)

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
ggsave(here::here("figures","figure6_E_coli.png"), width = 12, height = 7)

plot_vaccine_metric(data = chosen_results_to_plot,
                    metric_name_to_plot = c("prc_red_prop_prevCr"),
                    y_label = "Relative change due to vaccination (%)",
                    chosen_shape = chosen_shape,
                    chosen_palette = chosen_palette,
                    hrefs = c(0),
                    vrefs = c(30))
ggsave(here::here("figures",folder_name,paste0("prc_red_prop_prevCr_total.png")), width = 14, height = 5)


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

plot_vaccine_metric_both_bacteria(data1 = results_to_plot_Saureus,
                                  data2 = results_to_plot_Ecoli,
                                  metric_name_to_plot = c("prc_red_prop_prevCr"),
                                  y_label = "Relative change due to vaccination (%)",
                                  chosen_shapes = shapes,
                                  chosen_palette = palette_combined_named,
                                  hrefs = c(),
                                  vrefs = c(30),
                                  logscale = F)
ggsave(here::here("figures","prop_combined.png"), width = 14, height = 5)


#### Plot combined antibiotic impact ####

plot_antibiotic_metric_both_bacteria(data1 = results_to_plot_antibiotic_Saureus,
                                     data2 = results_to_plot_antibiotic_Ecoli,
                                     metric_name_to_plot = c("prc_red_inccumI", "prc_red_prop_prevCr"),
                                     chosen_shapes = shapes)


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


