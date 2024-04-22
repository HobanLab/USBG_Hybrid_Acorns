#####################
#     Libraries     #
#####################

library(diveRsity)
library(adegenet)
library(poppr)

###############################
#     Data Cleaning steps     #
###############################

#set working directoy
setwd("C:/Users/eschumacher/Desktop/UHA_Analysis")

#convert arp file to genepop file 
arp2gen("UHA_Final_Scores.arp")

#load in genepop file as a genind object
UHA_fs_genind <- read.genepop("UHA_Final_Scores.gen", ncode = 3)


#reduce genind file for individuals with greater than 25% missing data 
UHA_fs_genind_nomd <- missingno(UHA_fs_genind, type = "geno", 
                                cutoff = 0.25, quiet = FALSE, freq = FALSE)

#write out genind object as a genalex file
genind2genalex(UHA_fs_genind_nomd,
               "C:/Users/eschumacher/Desktop/UHA_Analysis/UHA_Final_Scores_nomd_genalex.csv")
