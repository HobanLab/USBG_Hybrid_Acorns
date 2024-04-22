#####################
#     Libraries     #
#####################

library(adegenet)
library(poppr)

###############################
#     Data Cleaning steps     #
###############################

#set working directoy
setwd("../Data_Files")

#load in genepop file as a genind object
UHA_genind <- read.genepop("Genotype_Files/2024_UHA_genepop.gen", ncode = 2)

#reduce genind file for individuals with greater than 25% missing data 
UHA_genind_nomd <- missingno(UHA_genind, type = "geno", 
                                cutoff = 0.25, quiet = FALSE, freq = FALSE)

#write out genind object as a genalex file
genind2genalex(UHA_genind_nomd,
               "CSV_Files/UHA_Final_Scores_clean.csv")
