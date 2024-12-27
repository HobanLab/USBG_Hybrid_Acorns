#####################
#     Libraries     #
#####################

library(tidyverse)

#######################
#     Load files      #
#######################

UHA_relate <- read.csv("./Data_Files/CSV_Files/UHA_spatial_distance_analysis.csv")

#####################################################
#     Reduced Individuals based on Relatedness      #
#####################################################

##reduce individuals based on 25% relatedness
#now identify how many individuals have greater than 25% relatedness = half siblings
UHA_halfsib_df <- UHA_relate %>%
                      dplyr::filter(relate_coef > 0.25)

#save out list of halfsib names 
UHA_halfsib_names <- unique(UHA_halfsib_df$Ind1)

#then use this to create a document which has all the unique individual numbers for every highly related individuals
butternut_halfsib_names_cleanfront <- gsub("^.*\\.","", butternut_halfsib_names)

butternut_halfsib_names_cleanback <- gsub("^.*\\_","", butternut_halfsib_names_cleanfront)
