#####################
#     Libraries     #
#####################

library(tidyverse)
library(geosphere)

######################
#     Load files     #
######################

#set working directory
setwd("../..")

#load in the tissue database, remove offspring which have no coordinates
UHA_db <- read.csv("Data_Files/CSV_Files/ARCHIVED_USBG_Hybrid_Acorn_Tissue_Database.csv")

#read in cleaned data file
UHA_score_clean_df <- read.csv("Data_Files/CSV_Files/UHA_score_clean_df.csv")  

#load in parentage results 
parentage_results <- read.csv("Data_Files/CSV_Files/UHA_full_parentage.csv")

###################################
#     Reorganize data frames      #
###################################
#remove individuals with too much missing data
UHA_clean_df <- UHA_db[UHA_db$Tissue_ID %in% UHA_score_clean_df$Tissue_ID,]

#remove offspring from the data frame
UHA_par_df <- UHA_clean_df %>% 
                filter(!is.na(Longitude) | !is.na(Latitude))
