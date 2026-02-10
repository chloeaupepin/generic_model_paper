
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

