###############################

start <- Sys.time()
source("read_header.R")

##### (1) functions ###########


#delete negative numbers
h <- function(x) {
  ifelse(x < 0, NA,  x)
}
#delete in y if x is over 10
k <- function(x, y) {
  ifelse(x > 10, NA, y)
}


##### (2) identify all objects in folder ###########
#
getwd()
folders <- list.dirs(
  path = paste0(getwd()), recursive = T) 

p <- paste
pprint <- function(...) print(paste(...))

#################
list_files <- function(folder) {
  print(p("looking in ", folder))
  
  files <- list.files(folder, pattern = ".*ASC", full.names = TRUE)
  
  remove_bashlash <- function(p) { #make path readable for R and windows
    gsub("\\\\", "/", p)
  }
  files <- lapply(files, remove_bashlash)
  
  return(files)
}


files <- list_files(folders)

##### (3) Read all files and put them in one df


read_file <- function(fn) {
  pprint("reading", fn)
  
  temp <- readLines(fn)
  temp2 <- head(temp[-(1:19)], - 4)
  temp2 <- head(temp2[-(2:3)],length(temp2))
  writeLines(temp2, con = "temporary.ASC")
  rawdata <- read.table(text = gsub(",", "\t", readLines("temporary.ASC"))
                        , sep = "\t", dec = ".", header = TRUE)
  # delete negative numbers
  rawdata[, ] <- lapply(rawdata[, ], h)
  # 
####### (4) only keep column Isotope & Intensity.AVG ######

    # use intereference corrected intensities
  rawdata_reduced_notime_int <- select(rawdata, "Isotope", "IF.cor.AVG")
 
  
  
  nm <- sub(".*/", "", fn, perl = TRUE)
  nm <- sub(".ASC", "", nm)
  print("-------------------------------------------------------")
  print(names(rawdata_reduced_notime_int))
  
  
  
  
  names(rawdata_reduced_notime_int)[names(rawdata_reduced_notime_int) == "IF.cor.AVG"] <- nm
  ## add time stamp ####
  asc_to_dat = sub(".ASC", ".dat", fn)
  time_stamp = extract_timestamp(asc_to_dat)
  times = list("time", time_stamp)

  rawdata_reduced_int <- rbind(times, rawdata_reduced_notime_int)
  
  return(rawdata_reduced_int)

}
################################################################################

transpose_list <- function(l) {
  # transpose data with first column as row names
  row_names = l$Isotope
  row_names = sub("\\(LR\\)", "", row_names)
  
  id <- colnames(l)[[2]]

  
  transposed <- as.data.frame(t(l[, -1]))
  colnames(transposed) <- row_names
  
  transposed <- cbind(id = id, transposed)
  
  return(transposed)
}

try_read_transpose <- function(fn) {
  ret <- tryCatch(
    {transpose_list(read_file(fn))},
    error = function(e) {
      return(data.frame())
    }
  )
}

transposed_reduced <- lapply(files, try_read_transpose)

results_int_ms <- rbindlist(transposed_reduced, fill = TRUE)


colnames(results_int_ms)[3:ncol(results_int_ms)] <- gsub(".HR.", "", colnames(results_int_ms)[3:ncol(results_int_ms)])


#######################
end <- Sys.time()
time_taken <- end - start
print(time_taken)

