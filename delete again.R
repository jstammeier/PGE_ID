library(readxl)
library(magrittr)
library(dplyr)
library(data.table)
library(stringr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # rstudio only

########################## Variables and variable dataframes #############################

## Weights
weights_df = as.data.frame(read_excel("Einwagen.xlsx"))
weights_df$date = as.Date(weights_df$date, format = "%Y-%m-%d")

## Calculate dilution factors

weights_df = subset(weights_df, date < "2023-09-29")

weights_df = mutate(weights_df
                    ,
                    dilf_spike = EW_gesamt / EW_spike
                    ,
                    dilf_sample = EW_gesamt / EW_sample)

source("All_data_int.R")

# average intensities for each repeated measurement

# results_int_ms$avg_spike_abd_ref = colMeans(
#   results_int_ms[grep("^Spike_\\d"
#                       , results_int_ms$id), "abd_ref"])

########################## Constants #############################

## concentrations element in spike

spike_version_chr = "stock-solution_1 2023-09-01; working-solution_1; calibrated 2023-09-15"
conc_Pt_spike	 =	11.43
conc_Pd_spike  =	8.56
conc_Re_spike	 =  7.031
conc_Ir_spike	 =  1.59


########################## Calculations #############################
###### isotope ratios

#### Pt
# Spiked isotope    Pt-194
# Reference isotope Pt-195

#### Pd
# Spiked isotope    Pd-106
# Reference isotope Pd-105

#### Re
# Spiked isotope    Re-185
# Reference isotope Re-187

#### Ir
# Spiked isotope    Pt-191
# Reference isotope Pt-193

elements_ID_df = data.frame(
  element = c("Pt",	"Pd",	"Re",	"Ir"),
  ref_isotope	 = c("Pt195",	"Pd105",	"Re187",	"Ir193"),
  spike_isotope	 = c("Pt194",	"Pd106",	"Re185",	"Ir191")
)

try(drop(results_conc_ID))
results_conc_ID = select(results_int_ms, id)

for (i in elements_ID_df$element) {
  print(paste(i))
  #set respective isotopes for active element(i)
  ref_isotope = elements_ID_df[elements_ID_df$element == i, "ref_isotope"]
  spike_isotope = elements_ID_df[elements_ID_df$element == i, "spike_isotope"]
  
  ## calculate abundances of respecive isotopes
  # calculate the sum of all intensities per element(i)
  results_int_ms = results_int_ms %>%
    rowwise() %>%
    mutate(
      sum_of_int = sum( c_across( starts_with(
        as.character(paste(i))
      )
      ),
      na.rm = T), 
      # calculate abundance of respective isotopes
      abd_ref = c_across(all_of(ref_isotope)) / sum_of_int ,
      abd_spike = c_across(all_of(spike_isotope)) / sum_of_int ,
      ratio = c_across(all_of(spike_isotope)) / c_across(all_of(ref_isotope))
    ) %>% ungroup()
  print("###########################################################")
  print(colnames(results_int_ms))
  print(results_int_ms$abd_ref)
  
  # average of the natural abundance, averaged from all measurements denoted "calib" followed by a number,
  # i.e. the "normal" calibration standard
  results_int_ms$avg_natural_abd_ref = colMeans(results_int_ms[grep("^Calib_\\d", results_int_ms$id), "abd_ref"])
  results_int_ms$avg_natural_abd_spike = colMeans(results_int_ms[grep("^Calib_\\d", results_int_ms$id), "abd_spike"])
  
  # average of the isotope abundance in the spike, averaged from all measurements denoted "spike" followed by a number
  results_int_ms$avg_spike_abd_ref = colMeans(results_int_ms[grep("^Spike_\\d", results_int_ms$id), "abd_ref"])
  results_int_ms$avg_spike_abd_spike = colMeans(results_int_ms[grep("^Spike_\\d", results_int_ms$id), "abd_spike"])
  # FUCK sollte technisch gesehen hirarchisch gemittelt werden...
  
  results_int_ms = results_int_ms %>%
    rowwise() %>%
    mutate(conc = (avg_spike_abd_spike - (ratio * avg_spike_abd_ref)) / (( ratio * avg_natural_abd_ref) - avg_natural_abd_spike )
    )
  
  results_conc_ID[[i]] = results_int_ms$conc
  # As-RbBs
  # Rb*Bx-Ax
  
  names(results_int_ms)[names(results_int_ms) == 'sum_of_int'] <- paste0('sum_of_int_', i)
  names(results_int_ms)[names(results_int_ms) == 'abd_ref'] <- paste0('abd_ref_', ref_isotope)
  names(results_int_ms)[names(results_int_ms) == 'abd_spike'] <- paste0('abd_spike_', spike_isotope)
  names(results_int_ms)[names(results_int_ms) == 'ratio'] <- paste0(ref_isotope, "/", spike_isotope)
  names(results_int_ms)[names(results_int_ms) == 'avg_natural_abd_ref'] <- paste0('avg_natural_abd_', ref_isotope)
  names(results_int_ms)[names(results_int_ms) == 'avg_natural_abd_spike'] <- paste0('avg_natural_abd_', spike_isotope)
  names(results_int_ms)[names(results_int_ms) == 'avg_spike_abd_ref'] <- paste0('avg_spike_abd_', ref_isotope)
  names(results_int_ms)[names(results_int_ms) == 'avg_spike_abd_spike'] <- paste0('avg_spike_abd_spike_', spike_isotope)
  
  
  
  
  
}  

results_int_ms = results_int_ms[, -which(names(results_int_ms) %in% c(  "Ru100"
                                                                        , "Ru101"
                                                                        , "Ru102"
                                                                        , "Ru99"
                                                                        , "Ru98"
                                                                        , "Ru96"
                                                                        , "Ru104"
                                                                        , "In115"
                                                                        , "Os189"
                                                                        , "Mo95"
                                                                        , "Sn118"
                                                                        , "Zr90"
                                                                        , "Hg200"
                                                                        , "Cd111"
))]

# test = results_conc_ID[grep("^Calib_Spike_\\d", results_int_ms$id),]

View(test)
View(results_conc_ID)

