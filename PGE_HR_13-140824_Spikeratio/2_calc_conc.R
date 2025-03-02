library(readxl)
library(magrittr)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(anytime)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # rstudio only

######################### Constants #############################

## concentrations element in spike

spike_version_chr = "stock-solution_1 2023-09-01; working-solution_1; calibrated 2023-09-15"

########################## Variables and variable data frames #############################

## Weights
weights_df = as.data.frame(read_excel("./Einwagen.xlsx"))
weights_df$date = as.Date(weights_df$date, format = "%Y-%m-%d")

## subset targeted date for concentration calculation
calc_conc_start_str <- "2024-06-18"
calc_conc_end_str <- "2024-08-26"

calc_conc_start_str <- as.Date(reversecalib_start_str, format = "%Y-%m-%d")
calc_conc_end_str <- as.Date(reversecalib_end_str, format = "%Y-%m-%d")

weights_df <- filter(weights_df, date >= calc_conc_start_str & date <= calc_conc_end_str)



##### read raw data as intensity
source("All_data_int_IFcorrint_header.R")

## Change Epoch datetime to Posix formatting
results_int_ms$time <- anytime(results_int_ms$time ,  tz = "CET") 

# Change datatable to data frame

results_int_ms %>%
  data.frame() -> results_int_ms

results_int_ms <- results_int_ms[, (names(results_int_ms) %in% columns_to_keep)]

# average intensities for each repeated measurement. Repeated measurements MUST be denoted "-1234"

unify_names_f = function(x){

  x = gsub("\\-\\d+$","", x)

  return(x)
}

# results_int_ms$id = unify_names_f(results_int_ms$id)

# Use column names to select columns (excluding the first two) to convert everything to numeric
cols_to_convert <- names(results_int_ms)[4:ncol(results_int_ms)]

for (col in cols_to_convert) {
  set(results_int_ms, j = col, value = as.numeric(results_int_ms[[col]]))
} 

# Change datatable to data frame

results_int_ms %>%
  data.frame() -> results_int_ms

################## Data cleaning ##############


results_int_ms <- results_int_ms[, (names(results_int_ms) %in% columns_to_keep)]

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
  spike_isotope	 = c("Pt194",	"Pd106",	"Re185",	"Ir191"),
  concentration = conc_spike
)

###### Create a dataframe to store results for each sample and element ####
results_conc_individual <- data.frame()

###########################################################################

# Create a unique identifier by combining id, time, and filepath
results_int_ms <- results_int_ms %>%
  mutate(unique_id = paste(filepath, id, sep = "*"))



for (unique_sample_id in unique(results_int_ms$unique_id)) {
  # Extract the original sample_id if needed later
  sample_id <- unique(results_int_ms$id[results_int_ms$unique_id == unique_sample_id]) # Extract original sample_id
  
  
  # Only process samples containing "+Sp"
  if (!grepl("\\+Sp", sample_id)) {
    next  # Skip to the next sample if it doesn't contain "+Sp"
  }

  for (i in elements_ID_df$element) {
    print(paste("********************************* Begin Calculation *************************************"))
    print(paste("*****************************************************************************************"))
    print(paste("Calculating for sample:", sample_id, "element:", i))
    print(paste("*****************************************************************************************"))
    
    ref_isotope <- elements_ID_df[elements_ID_df$element == i, "ref_isotope"]
    spike_isotope <- elements_ID_df[elements_ID_df$element == i, "spike_isotope"]
    conc_element <- elements_ID_df[elements_ID_df$element == i, "concentration"]
    
    # Calculate abundances and ratios for all samples sample
    # this is done per element
    results_int_ms <- results_int_ms %>%
      rowwise() %>%
      mutate(
        sum_of_int = sum(c_across(starts_with(as.character(paste(i)))), na.rm = TRUE),
        abd_ref = c_across(all_of(ref_isotope)) / sum_of_int,
        abd_spike = c_across(all_of(spike_isotope)) / sum_of_int,
        ratio = c_across(all_of(ref_isotope)) / c_across(all_of(spike_isotope))
      ) %>% ungroup()
    
 
    # Filter data for the current sample.  
    results_int_ms_sample <- results_int_ms %>% filter(unique_id == unique_sample_id)
    
    # Extract the ratio of the blend
    ratio_blend <- results_int_ms_sample$ratio
 

    # Extract Natural Abundance (from matching sample without "+Sp" or QC)
    matching_sample_id <- gsub("\\+Sp", "", sample_id)  # Remove "+Sp"
    matching_sample <- results_int_ms %>% filter(id == matching_sample_id)

    
    if (nrow(matching_sample) > 0) {
      
      avg_natural_abd_ref <- colMeans(as.matrix(matching_sample[, "abd_ref", drop = FALSE]), na.rm = TRUE)
      avg_natural_abd_spike <- colMeans(as.matrix(matching_sample[, "abd_spike", drop = FALSE]), na.rm = TRUE)
      print("Natural abundance was used from matching sample")

      } else {
      qc_samples <- results_int_ms %>% filter(grepl("Calib_1ppb", id))
      print(paste( "Using Calib_1ppb reference material instead:", qc_samples$id))
      if (nrow(qc_samples) > 0) {
        avg_natural_abd_ref <- colMeans(as.matrix(qc_samples[, "abd_ref", drop = FALSE]), na.rm = TRUE)
        avg_natural_abd_spike <- colMeans(as.matrix(qc_samples[, "abd_spike", drop = FALSE]), na.rm = TRUE)

        } else {
        # If neither nor exists
        warning(
          paste(
            "No matching sample or QC samples found for",
            sample_id,
            ". Using NA for natural abundance."
          )
        )
        avg_natural_abd_ref <- NA
        avg_natural_abd_spike <- NA
      }
    }
    

    # Calculate Spike Abundance (from Blank+Sp)
    blank_sp_samples <- results_int_ms %>% filter(grepl("8_Blank\\+Sp", id))
    if (nrow(blank_sp_samples) > 0) {
      avg_spike_abd_ref <- colMeans(as.matrix(blank_sp_samples[, "abd_ref", drop = FALSE]), na.rm = TRUE)
      avg_spike_abd_spike <- colMeans(as.matrix(blank_sp_samples[, "abd_spike", drop = FALSE]), na.rm = TRUE)
      } else {
      warning("No Blank+Sp samples found. Using NA spike abundance. Please look for another Blank+Sp")
      avg_spike_abd_ref <- NA
      avg_spike_abd_spike <- NA
     }

    conc = conc_element * ((ratio_blend * avg_spike_abd_spike) - avg_spike_abd_ref) / (avg_natural_abd_ref - (ratio_blend * avg_natural_abd_spike))
     
    print("------------------------- concentration is calculated as -------------------------")
    print(paste(conc_element, "*", ratio_blend, " * ", avg_spike_abd_spike
                , " - ",avg_spike_abd_ref
                , " / ",avg_natural_abd_ref, " - ",ratio_blend
                , " * ",avg_natural_abd_spike ))
    print("---------------------------------------------------------------------------------")
    
    results_conc_individual <- bind_rows(results_conc_individual,
                                         data.frame(id = unique_sample_id,
                                                    element = i,
                                                    conc = conc))
    
    print(paste("################################## The End ####################################### "))
    
  }
}
    
# Reshape the data to have elements as columns
results_conc_individual <- results_conc_individual %>%
  pivot_wider(
    id_cols = id, # Columns that together form the unique identifier
    names_from = element,   # Column whose values become the new column names
    values_from = conc      # Column whose values populate the new columns
  )

results_conc_individual$id <- str_split_i(results_conc_individual$id, "\\*", 2)

results_conc_individual <- merge(results_conc_individual, weights_df, by.x = "id", by.y = "id", all.x = TRUE)



###########

# step 1: resolve every dilution individually

results_conc = select(results_conc_individual, id)

for(j in elements_ID_df$element){

    results_conc[j] = results_conc_individual[j] * results_conc_individual$EW_spike / results_conc_individual$EW_sample

}



combine_replicates_f = function(x){
  x = gsub("\\_\\d$","", x)
  
  return(x)
}

# step 2: now calculate the mean for a sample across all dilutions. Separate dilutions MUST be denoted "_1234"

results_conc$id = combine_replicates_f(results_conc$id)
xxx
results_conc %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(
    n = n(),  
    
    mean_Pt = mean(Pt, na.rm = TRUE),
    mean_Pd = mean(Pd, na.rm = TRUE),    
    mean_Re = mean(Re, na.rm = TRUE),
    mean_Ir = mean(Ir, na.rm = TRUE),
    
    sd_Pt = sd(Pt, na.rm = TRUE),
    sd_Pd = sd(Pd, na.rm = TRUE),     
    sd_Re = sd(Re, na.rm = TRUE),
    sd_Ir = sd(Ir, na.rm = TRUE),
    
  ) -> results_conc


View(results_conc)


