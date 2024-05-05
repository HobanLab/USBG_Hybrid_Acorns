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


#par_genotyoes <- read.csv("UHA_score_database_nomd_genotype - UHA_score_database_nomd_genotype.csv") #read in genotype file to use species information of maternal and paternal trees.
#tissue_info <- read.csv("TCB_Tissue_Database.csv")

#################################################
#     Data Cleaning Post Parentage Analysis     #
#################################################

###add parental species to the data frame 

#joining the genotypes with the parentage results 
#to create the full_parentage data set that will be have spatial data, 
#hybrid status, and half sibling status.
full_parentage <- left_join(off_df, par_results, by='Offspring_ID') 

#Adding in Species Information for Maternal and Paternal trees
keep_col_ID <- c("Offspring_ID","Mother_ID", "Candidate_father_ID", "Species")

# Narrowing the data again after the join of data sets
full_parentage <- full_parentage[keep_col_ID] 
full_parentage <- full_parentage %>%
                    rename('Mother_Species' = 'Species')  

# #create data frame with father species 
# paternal_assignment <-  par_df[par_df$Tissue_ID %in% 
#                                  par_results$Candidate_father_ID,]

#add paternal species 
full_parentage <- left_join(full_parentage, par_df, 
                             by=c('Candidate_father_ID' = 'Tissue_ID'))  

#now reduce columns 
keep_col_ID2 <- c("Offspring_ID","Mother_ID","Candidate_father_ID", "Mother_Species", "Species")

#reduce data frame by populated columns
full_parentage <- full_parentage[keep_col_ID2]
#rename columns 
full_parentage <- full_parentage %>%
                    rename('Candidate_Father_Species' = 'Species')

###add geographic information for both parents to the data frame 

# #Adding in Maternal Accession and Longitude/Latitude information from TCB Tissue Database
# ```{r}
full_parentage <- read_csv("UHA_full_parentage.csv")
archived_tissues <- read_csv("ARCHIVED_USBG_Hybrid_Acorn_Tissue_Database.csv")
full_parentage <- left_join(full_parentage, archived_tissues, by=c('Mother ID' = 'Tissue ID'))
keeps <- c("Offspring ID","Mother ID", "Candidate father ID", "Mother Species", "Candidate Father Species", "Mother Accession", "Longitude", "Latitude")
full_parentage <- full_parentage[keeps] # Narrowing the data again after the join of data sets
full_parentage <- full_parentage %>%
  rename('Maternal Longitude' = 'Longitude', 'Maternal Latitude' = 'Latitude')
# # This is giving me some issues when I join data sets, sometimes the columns are titled with .x at the end and this affects the keeps values. If you rerun this code the keeps values might need to be adjusted
# ```
# 
# Adding in Candidate Father Accession and Longitude/Latitude information from TCB Tissue Database
# ```{r}
full_parentage <- left_join(full_parentage, archived_tissues, by=c('Candidate father ID' = 'Tissue ID'))
keeps <- c("Offspring ID","Mother ID", "Candidate father ID", "Mother Species", "Candidate Father Species", "Mother Accession", "Maternal Longitude", "Maternal Latitude", "Accession Number", "Longitude", "Latitude")
full_parentage <- full_parentage[keeps] # Narrowing the data again after the join of data sets
full_parentage <- full_parentage %>%
  rename('Candidate Father Accession' = 'Accession Number', 'Candidate Father Longitude' = 'Longitude', 'Candidate Father Latitude' = 'Latitude')
# # This has the same issues, sometimes the columns are titled with .x at the end and this affects the keeps values. If you rerun this code the keeps value might need to be adjusted
# ```
# 
# 
# Assigning Half Siblings in the data set
# ```{r}
# Mikaely Evans code for creating a new column to assign half sibling status to all the offspring
full_parentage$'Half Siblings' <- NA  # Made three new columns for this analysis
full_parentage$M_Accession_Abrv <- NA  
full_parentage$F_Accession_Abrv <- NA
# Abbreviating the accession numbers for the maternal and paternal trees was necessary to assign half sibling status because it is only necessary to look at the first 6 characters in the accession number to know if the trees came from the same lineage.

full_parentage$M_Accession_Abrv <- substr(full_parentage$'Mother Accession', 0, 6)
# This line adds the Maternal accession abbreviation to the M_Accession_Abrv column

full_parentage$F_Accession_Abrv <- substr(full_parentage$'Candidate Father Accession', 0, 6)
# This line adds the paternal accession abbreviation to the F_Accession_Abrv column

full_parentage <- full_parentage %>%
  mutate('Half Siblings' = case_when(M_Accession_Abrv == F_Accession_Abrv ~ TRUE,
                                     M_Accession_Abrv != F_Accession_Abrv ~ FALSE))
# This chunk above uses mutate to change the 'Half Siblings' column to represent the cases when the maternal accession and paternal accession match, and when they don't. They are represented by short phrases that are easier for readers to understand when they are graphed below.
# 
# ```
# 
# 
# Creating Distance Matrix
# ```{r}c
## Emily Schumacher code for creating distance matrix
# Mikaely Evans edited to use for the full_parentage data set
##calculate distances
#create a column for distance between mom and dad 
full_parentage$distance_between_parents <- NA

#loop to calculate distance between parents
for(d in 1:nrow(full_parentage)){
  
  full_parentage$distance_between_parents[d] <- distm(full_parentage[d,7:8], 
                                                      full_parentage[d,10:11],
                                                      fun=distGeo)
  
}

#calculate mean distance between parents 
UHA_dist_matrix <- matrix(nrow = length(unique(full_parentage$'Mother ID')),
                          ncol = 1)
for(m in 1:length(unique(full_parentage$'Mother ID'))){
  
  UHA_dist_matrix[m,1] <- mean(full_parentage[full_parentage$'Mother ID' == unique(full_parentage$'Mother ID')[[m]],][,11])
  
}

#name matrix 
rownames(UHA_dist_matrix) <- unique(full_parentage$'Mother ID')
colnames(UHA_dist_matrix) <- "Mean_Dist_Parents"
```

Creating Hybrid Offspring Column
```{r}
full_parentage$Hybrid <- NA
full_parentage <- full_parentage %>%
  mutate(Hybrid = case_when(`Mother Species` == `Candidate Father Species` ~ FALSE,
                            `Mother Species` != `Candidate Father Species` ~ TRUE))
```


Writing full_parentage into a csv
```{r}
write.csv(full_parentage, "UHA_full_parentage.csv", row.names=FALSE)
```


