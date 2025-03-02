library(readxl)
library(magrittr)
library(dplyr)
library(data.table)
library(stringr)
library(anytime)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # rstudio only

########################## Constants #############################

## concentrations element in spike

spike_version_chr = "stock-solution_1 2023-09-01; working-solution_1; calibrated 2023-09-15"

########################## Variables and variable data frames #############################

## Weights
weights_df = as.data.frame(read_excel("./Einwagen.xlsx"))
weights_df$date = as.Date(weights_df$date, format = "%Y-%m-%d")

## subset targeted date for reverse calibration
reversecalib_start_str <- "2024-08-12"
reversecalib_end_str <- "2024-08-12"

reversecalib_start_date <- as.Date(reversecalib_start_str, format = "%Y-%m-%d")
reversecalib_end_date <- as.Date(reversecalib_end_str, format = "%Y-%m-%d")

weights_df <- filter(weights_df, date >= reversecalib_start_date & date <= reversecalib_end_date)


######################### read raw data as intensity & convert columns if necessary ######

source("All_data_int_IFcorrint_header.R")

## Change Epoch datetime to Posix formatting
results_int_ms$time <- anytime(results_int_ms$time ,  tz = "CET") 


# average intensities for each repeated measurement. Repeated measurements MUST be denoted "-1234"

unify_names_f = function(x){
  
  x = gsub("\\-\\d+","", x)

  return(x)
}

results_int_ms$id = unify_names_f(results_int_ms$id)

# Use column names to select columns (excluding the first two) to convert everything to numeric
cols_to_convert <- names(results_int_ms)[4:ncol(results_int_ms)]

for (col in cols_to_convert) {
  set(results_int_ms, j = col, value = as.numeric(results_int_ms[[col]]))
} 

# Change datatable to data frame

results_int_ms %>%
  # group_by(id) %>%
  # summarise_all(mean) %>% #uncomment if repeated measurements need to be averaged
  data.frame() -> results_int_ms

################## Data cleaning ##############

# only keep relevant columns
columns_to_keep <- c("id", "time"
                     , "filepath"
                     , "Pt194", "Pt195"
                     , "Pd105", "Pd106"
                     , "Re185", "Re187"
                     , "Ir191", "Ir193"
)

results_int_ms <- results_int_ms[, (names(results_int_ms) %in% columns_to_keep)]

# # Find rows where the 'id' column contains "QC" 
# # deleted because they have weird abundances
# rows_to_remove <- grepl("QC", results_int_ms$id)
# 
# # Remove the identified rows using negative indexing
# results_int_ms <- results_int_ms[!rows_to_remove, ]


########################## Calculations #############################
###### constants and definitions

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

# This is the composition of the spike solution.
elements_ID_df = data.frame(
  element = c("Pt", "Pd", "Re", "Ir"),
  ref_isotope	 = c("Pt194", "Pd105", "Re187", "Ir193"),
  spike_isotope	 = c("Pt195", "Pd106", "Re185", "Ir191"),
  concentration = c(10, 10, 10, 10)
)


try(drop(results_conc_ID))
results_conc_ID = select(results_int_ms, id)

########## abundances and ratios #####

for (i in elements_ID_df$element) {
  print(paste(i))
  #set respective isotopes for active element(i)
  ref_isotope = elements_ID_df[elements_ID_df$element == i, "ref_isotope"]
  spike_isotope = elements_ID_df[elements_ID_df$element == i, "spike_isotope"]
  conc_element = elements_ID_df[elements_ID_df$element == i, "concentration"]
   
  print(ref_isotope)
  print(spike_isotope)
  print(conc_element)
  
  
  ## calculate abundances of respective isotopes
  # calculate the sum of all intensities per element(i)
  results_int_ms = results_int_ms %>%
    rowwise() %>%
    dplyr::mutate(
      sum_of_int = sum( c_across( starts_with(
                                              as.character(paste(i))
                                              )
                                  ),
                        na.rm = T),
      # calculate abundance of respective isotopes
      abd_ref = c_across( all_of(ref_isotope)) / sum_of_int ,
      abd_spike = c_across(all_of(spike_isotope)) / sum_of_int ,
      ratio =  c_across(all_of(ref_isotope)) / c_across(all_of(spike_isotope))
    ) %>% ungroup()

  print("###########################################################")
  print(colnames(results_int_ms))
  print(results_int_ms$abd_ref)

  # average of the natural abundance, averaged from all measurements denoted "Labking L1",
  # i.e. the "normal" reference standard that is here used as assay standard, i.e. for calibration: Calib
  results_int_ms$avg_natural_abd_ref = colMeans(results_int_ms[grep("^L1", results_int_ms$id), "abd_ref"])
  results_int_ms$avg_natural_abd_spike = colMeans(results_int_ms[grep("^L1", results_int_ms$id), "abd_spike"])

  # average of the isotope abundance in the spike, averaged from all measurements denoted "Blank+Sp" 
  results_int_ms$avg_spike_abd_ref = colMeans(results_int_ms[grep("Blank\\+Sp", results_int_ms$id), "abd_ref"])
  results_int_ms$avg_spike_abd_spike = colMeans(results_int_ms[grep("Blank\\+Sp", results_int_ms$id), "abd_spike"])
  # FUCK sollte technisch gesehen hirarchisch gemittelt werden... ist aber der Fall wenn es nur ein Run ist

  results_int_ms = results_int_ms %>%
    rowwise() %>%
    mutate(
      conc = conc_element *  (avg_natural_abd_ref - ( ratio * avg_natural_abd_spike )) / ((ratio * avg_spike_abd_spike) - avg_spike_abd_ref)
    )

  results_conc_ID[[i]] = results_int_ms$conc


  names(results_int_ms)[names(results_int_ms) == 'sum_of_int'] <- paste0('sum_of_int_', i)
  names(results_int_ms)[names(results_int_ms) == 'abd_ref'] <- paste0('abd_ref_', ref_isotope)
  names(results_int_ms)[names(results_int_ms) == 'abd_spike'] <- paste0('abd_spike_', spike_isotope)
  names(results_int_ms)[names(results_int_ms) == 'ratio'] <- paste0(ref_isotope, "/", spike_isotope)
  names(results_int_ms)[names(results_int_ms) == 'avg_natural_abd_ref'] <- paste0('avg_natural_abd_', ref_isotope)
  names(results_int_ms)[names(results_int_ms) == 'avg_natural_abd_spike'] <- paste0('avg_natural_abd_', spike_isotope)
  names(results_int_ms)[names(results_int_ms) == 'avg_spike_abd_ref'] <- paste0('avg_spike_abd_', ref_isotope)
  names(results_int_ms)[names(results_int_ms) == 'avg_spike_abd_spike'] <- paste0('avg_spike_abd_spike_', spike_isotope)


}  
  
######

results_conc_ID = merge(results_conc_ID, weights_df, by.x = "id", by.y = "id", incomparables = "NAN") 

# step 1: resolve every dilution individually

results_reversecalib_conc = select(results_conc_ID, id)

for(j in elements_ID_df$element){
  
  results_reversecalib_conc[j] = results_conc_ID[j] * results_conc_ID$EW_sample / results_conc_ID$EW_spike
  
}


combine_replicates_f = function(x){
  x = gsub("\\_\\d$","", x)
  
  return(x)
}

# step 2: now calculate the mean for a sample across all dilutions. Separate dilutions MUST be denoted "_1234"

results_reversecalib_conc$id = combine_replicates_f(results_reversecalib_conc$id)

results_reversecalib_conc %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(
    
    n = n(),  # Add the count of observations within each group
    
    mean_Pt = mean(Pt, na.rm = TRUE),
    mean_Pd = mean(Pd, na.rm = TRUE),    # Calculate mean for Pd, handling NAs
    mean_Re = mean(Re, na.rm = TRUE),
    mean_Ir = mean(Ir, na.rm = TRUE),
    
    sd_Pt = sd(Pt, na.rm = TRUE),
    sd_Pd = sd(Pd, na.rm = TRUE),     # Calculate sd for Pd, handling NAs
    sd_Re = sd(Re, na.rm = TRUE),
    sd_Ir = sd(Ir, na.rm = TRUE),
    
  ) -> results_reversecalib_conc


target_ids <- c("Calib+Sp1", "Calib+Sp2") 

conc_spike <- results_reversecalib_conc %>%
  filter(id == target_id) %>%
  select(mean_Pt, mean_Pd, mean_Re, mean_Ir) %>%
  unlist()

#Alternative using dplyr's `filter` with `%in%` and `group_by`

conc_spike <- results_reversecalib_conc %>%
  filter(id %in% target_ids) %>%
  summarise(
    mean_Pt = mean(mean_Pt, na.rm = TRUE),
    mean_Pd = mean(mean_Pd, na.rm = TRUE),
    mean_Re = mean(mean_Re, na.rm = TRUE),
    mean_Ir = mean(mean_Ir, na.rm = TRUE)
  ) %>%
  unlist()  


print(conc_spike)

View(results_reversecalib_conc)

# source("2_calc_conc.R")
