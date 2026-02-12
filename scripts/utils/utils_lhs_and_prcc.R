
# Transform lhs_matrix according to parameter distributions
transform_lhs_matrix <- function(lhs_matrix, param_names, param_distributions){
  lhs_transformed_matrix <- lhs_matrix
  for (i in seq_along(param_names)){
    quantile_function <- param_distributions[[param_names[i]]]
    lhs_transformed_matrix[, i] <- quantile_function(lhs_matrix[, i])
  }
  
  return(lhs_transformed_matrix)
}

# Generate a Latin Hypercube sampling associated to distribution given
generate_lhs <- function(number_of_samples, param_distributions, params_fixed_values){
  # list of parameter names
  param_names = names(param_distributions)
  
  # number of interesting parameters
  number_of_params = length(param_names)
  
  lhs_matrix <- randomLHS(number_of_samples, number_of_params)
  
  lhs_transformed_matrix <- transform_lhs_matrix(lhs_matrix, param_names, param_distributions)
  lhs_transformed_df <- as.data.frame(lhs_transformed_matrix)
  colnames(lhs_transformed_df) = param_names
  
  lhs_transformed_df <- merge(lhs_transformed_df, params_fixed_values)
  
  return(lhs_transformed_df)
}

#Compute prcc for one output
compute_prcc_epi <- function(input_df, output_vector){
  full_df <- bind_cols(input_df, tibble(output_name = output_vector))
  
  res <- epiR::epi.prcc(full_df, sided.test = 2, conf.level = 0.95)
  
  return(res)
}

#Comput prcc for multiple outputs
compute_prcc_epi_multi_outputs <- function(input_df, outputs,outputs_order){
  
  epi_prcc_results <- purrr::map(outputs, ~compute_prcc_epi(input_df, .x))
  
  #reformulation of results
  epi_prcc_df <- bind_rows(purrr::map(names(epi_prcc_results), function(output_name){
    prcc_output <- epi_prcc_results[[output_name]]
    tibble(
      output = output_name,
      param  = prcc_output$var,
      prcc   = prcc_output$est,
      lower  = prcc_output$lower,
      upper  = prcc_output$upper,
      pvalue = prcc_output$p.value)
  }) ) %>%
    mutate(output = factor(output, levels = outputs_order))
  
  epi_prcc_df$signif <- ifelse(epi_prcc_df$pvalue < 0.05, "Significant", "Non significant")
  return(epi_prcc_df)
}

output_labeller = c("prc_red_inccumI" = "Relative reduction of \ncumulative incidence \nof all infections",
                    "prc_red_inccumIs" = "Relative reduction of \ncumulative incidence \nof sensitive infections",
                    "prc_red_inccumIr" = "Relative reduction of \ncumulative incidence \nof resistant infections",
                    "inccumI" = "Cumulative incidence of \nall infections",
                    "inccumIs" = "Cumulative incidence of \nsensitive infections",
                    "inccumIr" = "Cumulative incidence of \nresistant infections",
                    "prop_inccumIr" = "Proportion of cumulative \nincidence of infections \ndue to resistant infections",
                    "prc_red_prop_inccumIr" = "Relative change of \nproportion of cumulative \nincidence of infections \ndue to resistant infections",
                    "prc_red_prevC" = "Relative reduction of \nprevalence of all colonizations",
                    "prc_red_prevCs" = "Relative reduction of \nprevalence of sensitive colonization",
                    "prc_red_prevCr" = "Relative reduction of \nprevalence of resistant colonization",
                    "prevC" = "Prevalence of \n all colonizations",
                    "prevCs" = "Prevalence of \nsensitive colonization",
                    "prevCr" = "Prevalence of \nresistant colonization",
                    "prop_prevCr" = "Proportion of prevalence \nof colonization \ndue to resistant colonization",
                    "prc_red_prop_prevCr" = "Relative change of \nproportion of prevalence \nof colonization \ndue to resistant colonization",
                    "Coexistence condition" = "Coexistence condition"

                    )

param_labeller = c("betaC" = bquote("\u03B2"[C]),
                   "betaI" = bquote("\u03B2"[I]),
                   "f"     = "f",
                   "as"    = "a",
                   "time_until_recovery_without_ATB_s" = bquote("\u03C4"[rec]),
                   "dps"   = "d",
                   "thetasr" = bquote("\u03B8"),
                   "prob_bystander_exposure"             = bquote(p[by]),
                   "time_until_decolo_by_bystander_ATB"  = bquote("\u03C4"["dec,by"]),
                   "prob_minority_strain_when_colonised" = bquote(p["min,C"]),
                   "prob_specific_exposure"              = bquote(p["sp,s"]),
                   "time_until_decolo_by_specific_ATB"   = bquote("\u03C4"["dec,sp,s"]),
                   "prob_minority_strain_when_infected" = bquote(p["min,I"]),
                   "prob_specific_exposure_r"              = bquote(p["sp,r"]),
                   "time_until_decolo_by_specific_ATB_r"   = bquote("\u03C4"["dec,sp,r"]),
                   "Vperc" = bquote(V["perc"]),
                   "vftcs" = bquote(v[a]),
                   "vfds"= bquote(v[d]),
                   "vfis" = bquote(v[i]),
                   "vfrs" = bquote(v[r]))

param_category_colors = c("bacteria" = "#E19BE1", "bystander antibiotic exposure" = "#FFCC31", "specific antibiotic exposure" = "#FF8000", "vaccine" = "#9BB4E9")

prepare_for_prcc_plot <- function(epi_prcc_df){
  epi_prcc_df %>%
    mutate(param_category = case_when(
      param %in% c("betaC","betaI","f","as","time_until_recovery_without_ATB_s","dps", "thetasr") ~ "bacteria",
      param %in% c("prob_bystander_exposure" ,"time_until_decolo_by_bystander_ATB","prob_minority_strain_when_colonised") ~ "bystander antibiotic exposure",
      param %in% c("prob_specific_exposure","time_until_decolo_by_specific_ATB","prob_minority_strain_when_infected" ,
                   "prob_specific_exposure_r" ,"time_until_decolo_by_specific_ATB_r") ~ "specific antibiotic exposure",
      param %in% c("Vperc","vftcs", "vfis", "vfds", "vfrs") ~"vaccine"
    ),
    param = factor(param, levels = c("vfrs", "vfis","vfds","vftcs","Vperc",
                                     "time_until_decolo_by_specific_ATB_r","prob_specific_exposure_r" ,
                                     "prob_minority_strain_when_infected","time_until_decolo_by_specific_ATB", "prob_specific_exposure",
                                     "prob_minority_strain_when_colonised","time_until_decolo_by_bystander_ATB","prob_bystander_exposure" ,
                                     "time_until_recovery_without_ATB_s","thetasr","as","dps","f","betaI","betaC"))
    )
}

# plots with color on parameters
plot_prcc_epi_multi_color <- function(epi_prcc_df, column_nbr, chosen_output_labeller = output_labeller,
                                      chosen_param_category_colors = param_category_colors,chosen_param_labeller = param_labeller){
  
  epi_prcc_df <- prepare_for_prcc_plot(epi_prcc_df)
  
  p_epi_1 <- epi_prcc_df %>%
    filter(output %in% c("prc_red_inccumI","prc_red_inccumIs","prc_red_inccumIr","inccumI","inccumIs","inccumIr","prop_inccumIr", "prc_red_prop_inccumIr")) %>%
    ggplot(aes(x = param, y = prcc, fill = param_category)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    ylim(-1,1)+
    coord_flip() +
    facet_wrap(~output, ncol = column_nbr, 
               labeller = as_labeller(chosen_output_labeller))+
    scale_fill_manual(values = chosen_param_category_colors) +
    scale_x_discrete(labels=chosen_param_labeller)+
    labs(
      title = "Cumulative incidence of infections",
      x = "Parameters",#"Variables d'entrée",
      y = "PRCC coefficient", #"Coefficient PRCC",
      fill = "Parameter category:"
    ) +
    theme_bw()  + theme(legend.position = "bottom",
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 14),
                        legend.text = element_text(size = 14),
                        legend.title = element_text(size=14),
                        plot.title=element_text(size=16),
                        strip.text= element_text(size=14))
  
  p_epi_2 <- epi_prcc_df %>%
    filter(output %in% c("prc_red_prevC","prc_red_prevCs","prc_red_prevCr","prevC","prevCs","prevCr","prop_prevCr","prc_red_prop_prevCr")) %>%
    ggplot(aes(x = param, y = prcc, fill = param_category)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    ylim(-1,1)+
    coord_flip() +
    facet_wrap(~output, ncol = column_nbr, 
               labeller = as_labeller(chosen_output_labeller))+
    scale_fill_manual(values = chosen_param_category_colors) +
    scale_x_discrete(labels=chosen_param_labeller)+
    labs(
      title = "Prevalence of colonization",
      x = "Parameters",#"Variables d'entrée",
      y = "PRCC coefficient", #"Coefficient PRCC",
      fill = "Parameter category:"
    ) +
    theme_bw() + theme(legend.position = "bottom",
                       axis.text = element_text(size = 12),
                       axis.title = element_text(size = 14),
                       legend.text = element_text(size = 14),
                       legend.title = element_text(size=14),
                       plot.title=element_text(size=16),
                       strip.text= element_text(size=14))
  
  p_epi = p_epi_1 / p_epi_2 + plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") & theme(legend.position = "bottom")
  return(p_epi)
}

plot_prcc_epi_multi_color_only_inccum_or_prev <- function(epi_prcc_df, column_nbr,cols_choice = c("prc_red_inccumI","prc_red_inccumIs","prc_red_inccumIr","inccumI","inccumIs","inccumIr","prop_inccumIr","prc_red_prop_inccumIr") ,chosen_output_labeller = output_labeller,
                                                  chosen_param_category_colors = param_category_colors,chosen_param_labeller = param_labeller){
  
  epi_prcc_df <- prepare_for_prcc_plot(epi_prcc_df)
  
  p_epi_1 <- epi_prcc_df %>%
    filter(output %in% cols_choice) %>%
    ggplot(aes(x = param, y = prcc, fill = param_category)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    ylim(-1,1)+
    coord_flip() +
    facet_wrap(~output, ncol = column_nbr, 
               labeller = as_labeller(chosen_output_labeller))+
    scale_fill_manual(values = chosen_param_category_colors) +
    scale_x_discrete(labels=chosen_param_labeller)+
    labs(
      x = "Parameters",#"Variables d'entrée",
      y = "PRCC coefficient", #"Coefficient PRCC",
      fill = "Parameter category:"
    ) +
    theme_bw()  + theme(axis.text = element_text(size = 12),
                        axis.text.y = element_text(size=14),
                        axis.title = element_text(size = 14),
                        legend.text = element_text(size = 14),
                        legend.title = element_text(size=14),
                        plot.title=element_text(size=16),
                        strip.text= element_text(size=14))
  
  
  
  p_epi = p_epi_1 + plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") 
  return(p_epi)
}
