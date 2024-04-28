#####################
#     Libraries     #
#####################

library(adegenet)
library(poppr)
library(PopGenReport)

###########################
#     Load Data Files     #
###########################

#set working directoy
setwd("../Data_Files")

#load in genepop file as a genind object
UHA_genind <- read.genepop("Genotype_Files/2024_UHA_genepop.gen", ncode = 2)

#load CERVUS data files
UHA_CERVUS_prep <- read.csv("Parentage_Files/2024_04_22_CERVUS_Prep.csv")

###############################
#     Data Cleaning steps     #
###############################

#reduce genind file for individuals with greater than 25% missing data 
UHA_genind_nomd <- missingno(UHA_genind, type = "geno", 
                                cutoff = 0.25, quiet = FALSE, freq = FALSE)

#write out genind object as a genalex file
genind2genalex(UHA_genind_nomd,
               "CSV_Files/UHA_Final_Scores_clean.csv")

###########################
#     Score Analysis      #
###########################

#run null allele calculations over all genind objects
UHA_null_all <- null.all(UHA_genind_nomd)

#visualize 
null_all_df <- signif(data.frame(UHA_null_all$null.allele.freq$summary2),3)*100

#calculate linkage disequilbrium
UHA_ld <- pair.ia(UHA_genind_nomd, sample = 1000)
#convert to a data frame
UHA_ld_df <- data.frame(round(UHA_ld, digits = 2))
