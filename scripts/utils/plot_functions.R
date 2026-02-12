
# Personal theme for ggplot 
personal_theme <- theme_bw()+
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.subtitle = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# Personal colors
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

palette_combined_named = c(
  "0.3 - Bacteria 1" = palette_orange_dark[1],
  "0.6 - Bacteria 1" = palette_orange_dark[3],
  "0.9 - Bacteria 1" = palette_orange_dark[5],
  "0.3 - Bacteria 2" = palette_pink_purple[1],
  "0.6 - Bacteria 2" = palette_pink_purple[3],
  "0.9 - Bacteria 2" = palette_pink_purple[5]
)

shapes = c(21,23)

colors_and_shapes <- list(
  "S_aureus" = list(palette = palette_orange_dark_3, shape = shapes[1]),
  "E_coli" = list(palette = palette_pink_purple_3, shape = shapes[2])
)

# Plot histogram of equilibrium results

prepare_eq_results_for_plot <- function(eq_results){
  eq_results %>%
    mutate(total_colonized = (eq_Csnv + eq_Crnv)/100000,
           resistance_ratio = eq_Crnv / (eq_Csnv + eq_Crnv)) %>%
    pivot_longer(cols = c(total_colonized, resistance_ratio), names_to = "metric", values_to = "value") %>%
    mutate(metric_name = case_when(
      metric == "total_colonized" ~ "Proportion colonized",
      metric == "resistance_ratio" ~ "Proportion of resistant\namong colonized"
    )) %>%
    mutate(metric_name = factor(metric_name, levels = c("Proportion colonized", "Proportion of resistant\namong colonized")))
}

plot_histogram_eq_results <- function(eq_results, limits_prop_colonized, limits_resistance_ratio){
  prepare_eq_results_for_plot(eq_results) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30) +
    facet_wrap(~metric_name, scales = "free_x")+
    labs(x = "Values",
         y = "Count") +
    personal_theme+
    facetted_pos_scales(x = list(
      scale_x_continuous(limits = limits_prop_colonized),
      scale_x_continuous(limits = limits_resistance_ratio)
    ))
}

plot_histogram_eq_results_two_bacteria <- function(eq_results_1, eq_results_2){
  bind_rows(
    eq_results_1 %>% mutate(bacteria_name = "Bacteria 1"),
    eq_results_2 %>% mutate(bacteria_name = "Bacteria 2")
  ) %>%
    prepare_eq_results_for_plot() %>%
    ggplot(aes(x = value, fill = bacteria_name)) +
    geom_histogram(bins = 100) +
    facet_grid(bacteria_name~metric_name, scales = "free_x")+
    scale_fill_manual(values = c(palette_orange_dark[2], palette_pink_purple[2]))+
    labs(x = "Values",
         y = "Count") +
    personal_theme + theme(legend.position = "none")
}




# Lignes de référence horizontales
geom_hrefs <- function(y = c(0), linetype = "dashed", color = "grey") {
  lapply(y, function(yval) geom_hline(yintercept = yval, linetype = linetype, color = color))
}

# Ligne verticale de référence
geom_vrefs <- function(x = c(30), linetype = "dashed", color = "grey") {
  lapply(x, function(xval) geom_vline(xintercept = xval, linetype = linetype, color = color))
}

# Plot vaccine metric for each vaccine type 
plot_vaccine_metric <- function(data,
                                metric_name_to_plot,
                                population_to_plot = c("Total population","Vaccinated", "Non vaccinated"),
                                vaccine_effects_to_plot = c(0.3, 0.6, 0.9),
                                y_label,
                                chosen_shape,
                                chosen_palette,
                                ylim_values = NULL,
                                hrefs,
                                vrefs){
  p <- data %>%
    filter(metric_name %in% metric_name_to_plot, 
           population %in% population_to_plot) %>%
    filter(varying_param %in% vaccine_effects_to_plot) %>%
    ggplot(aes(x = Vperc*100)) +
    geom_ribbon(aes(ymin = q025, ymax = q975,
                    fill = varying_param,
                    group = varying_param),
                alpha = 0.2) +
    geom_hrefs(y = hrefs)+
    geom_vrefs(x = vrefs)+
    geom_point(aes(y = median, color = varying_param, fill = varying_param), 
               shape = chosen_shape, size = 2)+
    geom_line(aes(y = median, color = varying_param)) +
    {if (length(metric_name_to_plot)>1) facet_nested(population ~ name_renamed) else facet_nested(~ name_renamed)}+
    scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
    scale_color_manual(values = chosen_palette)+
    scale_fill_manual(values = chosen_palette)+
    labs(x = "Vaccine coverage (%)",
         y = y_label,
         color = "Vaccine efficacy",
         fill = "Vaccine efficacy")+
    personal_theme + 
    theme(legend.position = "bottom")+
    {if (!is.null(ylim_values)) ylim(ylim_values[1], ylim_values[2]) else NULL}
  
  return(p)
}


plot_vaccine_metric_both_bacteria <- function(data1,
                                              data2,
                                              metric_name_to_plot,
                                              population_to_plot = c("Total population","Vaccinated", "Non vaccinated"),
                                              vaccine_effects_to_plot = c(0.3, 0.6, 0.9),
                                              y_label,
                                              chosen_shapes,
                                              chosen_palette,
                                              ylim_values = NULL,
                                              hrefs,
                                              vrefs,
                                              logscale = F){
  
  data <- data1 %>%
    mutate(bacteria = "1") %>%
    bind_rows(data2 %>% mutate(bacteria = "2"))
  
  p <- data %>%
    filter(metric_name %in% metric_name_to_plot, 
           population %in% population_to_plot) %>%
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
    filter(varying_param %in% vaccine_effects_to_plot) %>%
    ggplot(aes(x = Vperc*100)) +
    geom_ribbon(aes(ymin = q025, ymax = q975,
                    fill = color_group,
                    group = color_group),
                alpha = 0.2) +
    geom_hrefs(y = hrefs)+
    geom_vrefs(x = vrefs)+
    geom_point(aes(y = median, color = color_group, fill = color_group, shape = color_group),size = 2)+
    geom_line(aes(y = median, color = color_group)) +
    {if (length(metric_name_to_plot)>1) facet_nested(population ~ name_renamed) else facet_nested(~ name_renamed)}+
    {if(logscale) scale_y_log10(minor_breaks = minor_breaks_log())}+
    scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
    scale_color_manual(values = chosen_palette)+
    scale_fill_manual(values = chosen_palette)+
    scale_shape_manual(values = c(chosen_shapes[1],chosen_shapes[1],chosen_shapes[1],chosen_shapes[2],chosen_shapes[2], chosen_shapes[2]))+
    labs(x = "Vaccine coverage (%)",
         y = y_label,
         color = "Vaccine efficacy",
         fill = "Vaccine efficacy")+
    personal_theme + 
    theme(legend.position = "bottom")+
    {if (!is.null(ylim_values)) ylim(ylim_values[1], ylim_values[2]) else NULL}
  
  p <- p + 
    guides(fill = guide_legend(nrow = 2, 
                               byrow = TRUE,
                               override.aes = list(shape = ifelse(grepl("Bacteria 1", levels(p$data$color_group)), 21, 23))), 
           color = guide_legend(nrow = 2, 
                                byrow = TRUE), 
           shape = "none")
  
  return(p)
}

prepare_for_antibiotic_metric_plot <- function(data, output_order){
  data %>% mutate(prc_red_prob_bystander_exposure = case_when(
    scenario_id == 1 ~ 10,
    scenario_id == 2 ~ 30,
    scenario_id == 3 ~ 50,
    scenario_id == 4 ~ 70,
    scenario_id == 5 ~ 90
  )) %>%
    mutate(metric_name = factor(metric_name, levels = output_order))
}


output_labeller = c("prc_red_inccumI" = "Cumulative incidence \nof all infections",
                    "prc_red_inccumIs" = "Cumulative incidence \nof sensitive infections",
                    "prc_red_inccumIr" = "Cumulative incidence \nof resistant infections",
                    "prc_red_prop_inccumIr" = "Proportion of cumulative \nincidence of infections \ndue to resistant infections",
                    "prc_red_prevC" = "Prevalence of \nall colonizations",
                    "prc_red_prevCs" = "Prevalence of \nsensitive colonization",
                    "prc_red_prevCr" = "Prevalence of \nresistant colonization",
                    "prc_red_prop_prevCr" = "Proportion of prevalence \nof colonization \ndue to resistant colonization"
                    
                    
)

plot_antibiotic_metric_both_bacteria <- function(data1,
                                   data2,
                                   metric_name_to_plot,
                                   chosen_shapes,
                                   facet_cols = NULL,
                                   chosen_output_labeller = output_labeller){
  
  data <- data1 %>%
    mutate(Bacteria = "1") %>%
    bind_rows(data2 %>% mutate(Bacteria = "2")) %>%
    prepare_for_antibiotic_metric_plot(., output_order = names(chosen_output_labeller))
  
  
  p <- data %>%
    filter(metric_name %in% metric_name_to_plot) %>%
    ggplot(aes(x=prc_red_prob_bystander_exposure)) + 
    geom_hline(yintercept = 0, linetype="dashed", color = "black")+
    geom_ribbon(aes(ymin = q025, ymax = q975, color = Bacteria, fill = Bacteria),
                alpha = 0.2) +
    geom_line(aes(y = median, color = Bacteria)) +
    geom_point(aes(y = median, color = Bacteria, fill = Bacteria, shape = Bacteria), size = 2)+
    {if (length(metric_name_to_plot)>1) facet_wrap(~ metric_name, ncol = facet_cols, 
                                                   labeller = as_labeller(chosen_output_labeller)) else NULL}+
    scale_x_continuous(breaks = c(10, 30, 50, 70, 90))+
    scale_color_manual(values = c(palette_orange_dark[3], palette_pink_purple[3]))+
    scale_fill_manual(values = c(palette_orange_dark[3], palette_pink_purple[3]))+
    scale_shape_manual(values = c(chosen_shapes[1],chosen_shapes[2]))+
    labs(
      x = "Reduction percentage of bystander antibiotic exposure (%)",
      y = "Relative change due to\nbystander antibiotic reduction (%)"
    ) +
    personal_theme +
    theme(legend.position = "bottom")
  
  return(p)
}

