library(here)
# Run the complete analysis for the two case studies 

file_name_S_aureus = "S_aureus_params_eu_with_J01CR.csv"
file_name_E_coli = "E_coli_params_primavera5.csv"
#### S. aureus ####

bacteria = "S_aureus"
folder_name = "S_aureus"
file_name  = file_name_S_aureus
transmission_by_infected = FALSE

source(here("scripts", "run_vaccination_impact_analysis.R"))
source(here("scripts","run_antibiotic_impact_analysis.R"))
source(here("scripts","vaccine_impact_on_R0.R"))

#### E. coli ####
bacteria = "E_coli"
folder_name = "E_coli"
file_name = file_name_E_coli
transmission_by_infected = TRUE # if transmission is allowed, it equals transmission by colonized individuals

source(here("scripts", "run_vaccination_impact_analysis.R"))
source(here("scripts","run_antibiotic_impact_analysis.R"))
source(here("scripts","vaccine_impact_on_R0.R"))

#### Supplementary analysis ####
source(here("scripts", "why_keep_f_constant.R"))
