###This code was updated on 5/5/24 by Emily Schumacher. This code was worked on by
##Mikaely Evans as well as Ash Hamilton to generate a data file, UHA_full_parentage
#which summarizes the results of parentage analysis and the remaining analyses of the project
#are run on. This data file combines information about the assigned father, 
#mother, distance between them, whether the offspring is a hybrid, 
#and how far apart the parents are, as well as if the parents are half-siblings. 

#####################
#     Libraries     #
#####################

library(tidyverse)
library(dplyr)
library(ggplot2)
library(geofi)
library(sf)
library(geosphere)

###########################
#     Load Data Files     #
###########################
#set working directory 
setwd("../..")

# parentage_results <- read_csv("UHA_11_14_23_parentage_analysis.csv")  # reading in parentage_results
#spatial_parentage_results <- read_csv("UHA_11_14_23_parentage_analysis.csv")  # reading in parentage results again for future addition of spatial data
#write in the parentage results summary 
par_results <- read.csv("Analysis/Parentage_Analysis/Initial_Run/Output_Files/UHA_parentage_sumary.csv")

#replace periods with underscones 
colnames(par_results) <-  gsub('\\.', '_', colnames(par_results))

#parentage df 
par_df <- read.csv("Data_Files/CSV_Files/UHA_score_clean_df.csv")
#remove first column
par_df <- par_df[,-1]

#rename offspring ID to tissue ID
colnames(off_df)[1] <- "Offspring_ID"

#limit to offspring 
off_df <- par_df[par_df$Parent_Offspring == "O",]

#load data file with geographic information 
UHA_database <- read.csv("Data_Files/CSV_Files/ARCHIVED_USBG_Hybrid_Acorn_Tissue_Database.csv")

#par_genotyoes <- read.csv("UHA_score_database_nomd_genotype - UHA_score_database_nomd_genotype.csv") #read in genotype file to use species information of maternal and paternal trees.
#tissue_info <- read.csv("TCB_Tissue_Database.csv")

#################################################
#     Data Cleaning Post Parentage Analysis     #
#################################################

###add parental species to the data frame 

#joining the genotypes with the parentage results 
#to create the full_parentage data set that will be have spatial data, 
#hybrid status, and half sibling status.
full_parentage <- left_join(off_df, par_results, by="Offspring_ID")

#Adding in Species Information for Maternal and Paternal trees
keep_col_ID <- c("Offspring_ID","Mother_ID", "Candidate_father_ID", "Species")

# Narrowing the data again after the join of data sets
full_parentage <- full_parentage[keep_col_ID] 
full_parentage <- full_parentage %>%
                    rename('Mother_Species' = 'Species')  

#add paternal species name
full_parentage <- left_join(full_parentage, par_df, 
                             by=c('Candidate_father_ID' = 'Tissue_ID'))  

#reduce by empty columns 
keep_col_ID2 <- c("Offspring_ID","Mother_ID","Candidate_father_ID", 
                  "Mother_Species", "Species")

#reduce data frame by populated columns
full_parentage <- full_parentage[keep_col_ID2]
#rename columns 
full_parentage <- full_parentage %>%
                    rename('Candidate_Father_Species' = 'Species')

###add geographic information for both parents to the data frame 

#Adding in Maternal Accession and Longitude/Latitude information from TCB Tissue Database
full_parentage <- left_join(full_parentage, UHA_database, 
                                    by=c('Mother_ID' = 'Tissue_ID'))

#rename accession column 
full_parentage <- full_parentage %>%
                    rename("Maternal_Accession" = "Accession.Number")

#remove added columns
keep_col_ID3 <- c("Offspring_ID", "Mother_ID", "Candidate_father_ID",
                  "Mother_Species", "Candidate_Father_Species",
                  "Maternal_Accession", "Longitude", "Latitude")
      
full_parentage <- full_parentage[keep_col_ID3] 

#rename longitude and latitude 
full_parentage <- full_parentage %>%
  rename(c('Maternal_Longitude' = 'Longitude', "Maternal_Latitude" = "Latitude"))

#Adding in Candidate Father Accession and Longitude/Latitude information from TCB Tissue Database
full_parentage <- left_join(full_parentage, 
                            UHA_database, by=c('Candidate_father_ID' = 'Tissue_ID'))

#column IDs to keep 
keep_col_ID4 <- c("Offspring_ID","Mother_ID", "Candidate_father_ID", 
                   "Mother_Species", "Candidate_Father_Species", "Maternal_Accession", 
                    "Maternal_Longitude", "Maternal_Latitude", "Accession.Number", 
                      "Longitude", "Latitude")

full_parentage <- full_parentage[keep_col_ID4] # Narrowing the data again after the join of data sets

full_parentage <- full_parentage %>%
  rename('Candidate_Father_Accession' = 'Accession.Number', 
         'Candidate_Father_Longitude' = 'Longitude', 
         'Candidate_Father_Latitude' = 'Latitude')

###Assigning Half Siblings in the data set
# Mikaely Evans code for creating a new column to assign half sibling status to all the offspring
full_parentage$Half_Sibs <- NA  # Made three new columns for this analysis
full_parentage$M_Accession_Abrv <- NA  
full_parentage$F_Accession_Abrv <- NA

#Abbreviating the accession numbers for the maternal and paternal trees was necessary to assign half sibling status because it is only necessary to look at the first 6 characters in the accession number to know if the trees came from the same lineage.
full_parentage$M_Accession_Abrv <- substr(full_parentage$Maternal_Accession, 0, 6)

#This line adds the Maternal accession abbreviation to the M_Accession_Abrv column
full_parentage$F_Accession_Abrv <- substr(full_parentage$Candidate_Father_Accession, 0, 6)

#This line adds the paternal accession abbreviation to the F_Accession_Abrv column
full_parentage <- full_parentage %>%
  mutate("Half_Sibs" = case_when(M_Accession_Abrv == F_Accession_Abrv ~ TRUE,
                                     M_Accession_Abrv != F_Accession_Abrv ~ FALSE))

# Creating Distance Matrix
# ```{r}c
## Emily Schumacher code for creating distance matrix
# Mikaely Evans edited to use for the full_parentage data set
##calculate distances
#create a column for distance between mom and dad 
full_parentage$dist_par <- NA

#loop to calculate distance between parents
for(d in 1:nrow(full_parentage)){
  
  full_parentage$dist_par[d] <- distm(full_parentage[d,7:8], 
                                            full_parentage[d,10:11],
                                                      fun=distGeo)
  
}

#calculate mean distance between parents 
UHA_dist_matrix <- matrix(nrow = length(unique(full_parentage$Mother_ID)),
                          ncol = 1)

for(m in 1:length(unique(full_parentage$Mother_ID))){
  
  UHA_dist_matrix[m,1] <- mean(full_parentage[full_parentage$Mother_ID == unique(full_parentage$Mother_ID)[[m]],][,11])
  
}

#name matrix 
rownames(UHA_dist_matrix) <- unique(full_parentage$Mother_ID)
colnames(UHA_dist_matrix) <- "Mean_Dist_Parents"

full_parentage$Hybrid_Status <- NA
full_parentage <- full_parentage %>%
  mutate(Hybrid_Status = case_when(Mother_Species == Candidate_Father_Species ~ FALSE,
                            Mother_Species != Candidate_Father_Species ~ TRUE))

#write to file 
write.csv(full_parentage, "Data_Files/CSV_Files/UHA_full_parentage.csv",
          row.names = FALSE)
