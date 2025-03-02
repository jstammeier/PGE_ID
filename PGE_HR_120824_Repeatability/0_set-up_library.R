setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # rstudio only
# found on https://statsandr.com/blog/an-efficient-way-to-install-and-load-r-packages/

# Package names
packages <- c(  "ggplot2"  , "tidyverse" , "ggQC"
              , "openxlsx" , "readxl"    , "magrittr"
              , "reshape2" , "see"       , "dplyr"
              , "patchwork", "plotly"    , "hrbrthemes"
              #til here necessary for xrf
              , "zoo"      , "hms"       , "ggplot2" #add necessary for import ms
              , "ggrepel"  , "anytime"   , "data.table"
              , "tools"    , "purrr"     , "stringr" 
              , "rlist"    , "reader"    , "knitr"    
              , "grid"     , "gridExtra" , "ggpubr"
              )


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, library, character.only = TRUE)


