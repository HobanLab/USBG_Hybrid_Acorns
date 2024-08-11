#This script was created to clean the occ records of all accessioned trees
#at the arb and limit to only the oak trees

#####################
#     Libraries     #
#####################

library(dplyr)
library(sjmisc)

#####################
#     Analysis      #
#####################

#set wd 
setwd("C:/Users/eschumacher/Documents/GitHub/USBG_Hybrid_Acorn")

#load in data file with all accessioned trees for the TMA 
TMA_all_trees <- read.csv("Data_Files/Geographic_Files/TMA_all_trees.csv")

#loop to ID which individuals are quercus records 
for(i in 1:nrow(TMA_all_trees)){

  TMA_all_trees$Quercus_tf[[i]] <- str_contains(TMA_all_trees$CalcFullName[[i]], pattern = "Quercus ",
                                         ignore.case = TRUE)
}

#spot check if the loop worked using this df 
TMA_all_trees_red <- TMA_all_trees[,c(1:6,40)]

#limit by quercus records 
TMA_quercus_occ <- TMA_all_trees[TMA_all_trees$Quercus_tf == TRUE,]

#write out csv
write.csv(TMA_quercus_occ,"Data_Files/Geographic_Files/TMA_quercus_occ.csv",
          row.names = FALSE)
