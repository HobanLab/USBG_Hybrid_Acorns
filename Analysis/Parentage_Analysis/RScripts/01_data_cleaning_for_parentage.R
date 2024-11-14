###This RScript processes genotyope files generated from geneious to
##clean for 

#####################
#     Libraries     #
#####################

library(adegenet)
library(poppr)
library(PopGenReport)
library(hierfstat)
library(tidyverse)

###########################
#     Load Data Files     #
###########################

#set working directoy
# setwd("../../..")

#load in genepop file as a genind object
UHA_genind <- read.genepop("Data_Files/Genotype_Files/UHA_raw_genotype.gen", ncode = 2)

#load score df
UHA_scores_df <- read.csv("Data_Files/CSV_Files/UHA_score_database.csv")

#load in UHA database
UHA_database <- read.csv("Data_Files/CSV_Files/UHA_database.csv")

###############################
#     Data Cleaning steps     #
###############################

#reduce genind file for individuals with greater than 25% missing data 
UHA_genind_nomd <- missingno(UHA_genind, type = "geno", 
                                cutoff = 0.25, quiet = FALSE, freq = FALSE)

#write out genind object as a genalex file
genind2genalex(UHA_genind_nomd,
               "Data_Files/Genotype_Files/UHA_genalex_clean.csv")

#limit by the cleaned individuals
UHA_scores_clean_df <- UHA_scores_df[UHA_scores_df[,1] %in% 
                                       rownames(UHA_genind_nomd@tab),]

#write out 
write.csv(UHA_scores_clean_df, "Data_Files/CSV_Files/UHA_score_clean_df.csv")


#######################################
#     Write out Summary Database      #
#######################################

#organized 
UHA_species <- as.data.frame(table(UHA_database$Species))

write.csv(UHA_species, "./Results/Preliminary_Genotyping_Analysis/UHA_species.csv")


###########################
#     Score Analysis      #
###########################

####Null alleles
#run null allele calculations over all genind objects
UHA_null_all <- null.all(UHA_genind_nomd)

#store in a data frame 
null_all_df <- signif(data.frame(UHA_null_all$null.allele.freq$summary2),3)*100

#write out null allele data frame
write.csv(null_all_df, "Results/Preliminary_Genotyping_Analysis/null_all_df.csv")

####Linkage Disequilibrium
#calculate linkage disequilbrium
UHA_ld <- pair.ia(UHA_genind_nomd, sample = 1000)

#convert to a data frame
UHA_ld_df <- data.frame(round(UHA_ld, digits = 2))

#write out ld df 
write.csv(UHA_ld_df, "Results/Preliminary_Genotyping_Analysis/ld_df.csv")

#########################
#     Pop Gen Stats     #
#########################

#create data table 
gendiv_sum_stats <- matrix(nrow = 3, ncol = 5)

#create poppr data file 
gendiv_poppr <- poppr(UHA_genind_nomd)

#number of individuals 
gendiv_sum_stats[1,1] <- gendiv_poppr$N[[1]]
gendiv_sum_stats[2,1] <- gendiv_poppr$N[[2]]
gendiv_sum_stats[3,1] <- gendiv_poppr$N[[3]]

#MLG
gendiv_sum_stats[1,2] <- gendiv_poppr$MLG[[1]]
gendiv_sum_stats[2,2] <- gendiv_poppr$MLG[[2]]
gendiv_sum_stats[3,2] <- gendiv_poppr$MLG[[3]]

#number of alleles 
gendiv_sum_stats[1,3] <- summary(UHA_genind_nomd)[[4]][[1]]
gendiv_sum_stats[2,3] <- summary(UHA_genind_nomd)[[4]][[2]]
gendiv_sum_stats[3,3] <- sum(nAll(UHA_genind_nomd))

#allrich
gendiv_sum_stats[1,4] <- mean(allelic.richness(UHA_genind_nomd)$Ar[,1])
gendiv_sum_stats[2,4] <- mean(allelic.richness(UHA_genind_nomd)$Ar[,2])
gendiv_sum_stats[3,4] <- mean(rowMeans(allelic.richness(UHA_genind_nomd)$Ar))

#expected heterozygosity
gendiv_sum_stats[1,5] <- gendiv_poppr$Hexp[[1]]
gendiv_sum_stats[2,5] <- gendiv_poppr$Hexp[[2]]
gendiv_sum_stats[3,5] <- gendiv_poppr$Hexp[[3]]

#name rows and columns 
rownames(gendiv_sum_stats) <- c("Offspring", "Parents", "Total")
colnames(gendiv_sum_stats) <- c("N", "MLG", "nAll", "Allelic Richness",
                                "Expected Heterozygosity")

#write out summary data file 
write.csv(gendiv_sum_stats, "Results/Preliminary_Genotyping_Analysis/gendiv_sum_stats.csv")


#########################################
#     Prep Data Files for Parentage     #
#########################################

##Create data frames with and without loci with 
#large number of null alleles

#first identify which alleles have > 10% null alleles
na_reorg_df <- as.data.frame(t(null_all_df))

#clean row names 
rownames(na_reorg_df) <- gsub("\\...*", "",rownames(na_reorg_df))

#list of loci with > 15% null alleles 
na_loci <- gsub("\\_.*","",rownames(na_reorg_df[which(na_reorg_df$`Observed frequency` > 15),]))


#generate input data frame for parentage - all loci
write.csv(UHA_scores_clean_df, "Analysis/Parentage_Analysis/CERVUS_Files/All_Loci/Input_Files/UHA_all_loci_clean_genotype_df.csv",
          row.names = FALSE)

#generate a data frame removing loci with high null allele %
UHA_red_loci_df <- UHA_scores_clean_df %>%
                    dplyr::select(-contains(na_loci))

#generate parentage input data without loci with high null allele %
write.csv(UHA_red_loci_df, "Analysis/Parentage_Analysis/CERVUS_Files/Red_Loci/Input_Files/UHA_red_loci_clean_genotype_df.csv",
          row.names = FALSE)

