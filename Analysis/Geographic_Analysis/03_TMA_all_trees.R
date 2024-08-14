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
setwd("../..")

#load in data file with all accessioned trees for the TMA 
TMA_all_trees <- openxlsx::read.xlsx("Data_Files/Geographic_Files/TMA_all_trees.xlsx")

#loop to ID which individuals are quercus records 
for(i in 1:nrow(TMA_all_trees)){

  TMA_all_trees$Quercus_tf[[i]] <- str_contains(TMA_all_trees$CalcFullName[[i]], pattern = "Quercus ",
                                         ignore.case = TRUE)
}

#limit by quercus records 
TMA_quercus_occ <- TMA_all_trees[TMA_all_trees$Quercus_tf == TRUE,]

#write out csv
openxlsx::write.xlsx(TMA_quercus_occ,"Data_Files/Geographic_Files/TMA_quercus_occ.xlsx")
