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

#load parentage results - replacing:
#par_results <- read.csv("Analysis/Parentage_Analysis/Initial_Run/Output_Files/UHA_parentage_sumary.csv")

#load parentage results - all loci = al
UHA_al_par <- read.csv("Analysis/Parentage_Analysis/All_Loci/Output_Files/all_loci_par_sum.csv",
                       row.names = NULL)

#remove periods from UHA data frame
colnames(UHA_al_par) <- gsub("\\.", "_", colnames(UHA_al_par))

#load parentage results - reduced loci = rl
UHA_rl_par <- read.csv("Analysis/Parentage_Analysis/Red_Loci/Output_Files/red_loci_par_sum.csv")

#remove periods from colnames
colnames(UHA_rl_par) <- gsub("\\.", "_", colnames(UHA_rl_par))

##create a list of the parentage results data frames
par_res_list <- list(UHA_al_par, UHA_rl_par)

#create a list of scenarios 
scen <- c("all_loci", "red_loci")

#loop to reduce the parental summary data frames for low confidence father assignments
for(sc in 1:length(scen)){
  
  temp_df <- par_res_list[[sc]]
  
  #first, join genotype files with parentage results 
  HCF_df  <- temp_df[((temp_df$Pair_LOD_score > 0) & (temp_df$Trio_LOD_score > 0)),]
  
  #write out data frame 
  write.csv(HCF_df, paste0("Data_Files/CSV_Files/", scen[[sc]], "_HCF_par_sum.csv"))
  
}
 
#load in all parentage summary data frames 
par_sum_list <- list.files(path = "Data_Files/CSV_Files/", pattern = "par_sum")

#reorder list 
par_sum_list <- list(par_sum_list[[2]], par_sum_list[[1]], 
                     par_sum_list[[4]], par_sum_list[[3]])

#read CSVs 
par_sum_df_list <- list()

#Loop over 
for(df in par_sum_list){
  
  par_sum_df_list[[df]] <- read.csv(paste0("Data_Files/CSV_Files/", df))
  
}

##load score data frames
#All loci
al_score_df <- read.csv("Analysis/Parentage_Analysis/All_Loci/Input_Files/UHA_all_loci_genotype_df.csv")

#Reduced loci
rl_score_df <- read.csv("Analysis/Parentage_Analysis/Red_Loci/Input_Files/UHA_red_loci_genotype_df.csv")

#combine into a list 
score_df_list <- list(al_score_df, rl_score_df)


###create offspring only data frames
#list of offspring data frames

off_df_list <- list()

for(o in 1:length(score_df_list)){
  
  #limit data frame to just offspring for each scenario 
  off_df_list[[o]] <- score_df_list[[o]][score_df_list[[o]]$Parent_Offspring == "O",]
  
  #rename "Tissue ID" to "Offspring_ID"
  colnames(off_df_list[[o]])[1] <- "Offspring_ID"
  
}

#load UHA database 
UHA_database <- read.csv("Data_Files/CSV_Files/ARCHIVED_USBG_Hybrid_Acorn_Tissue_Database.csv")

######### Create analysis data frame -------------------

#all scenarios 
full_scen <- c("all_loci", "HCF_all_loci", 
               "red_loci", "HCF_red_loci")
#loop over four scenarios
for(sc in 1:length(full_scen)){
  
  #store temporary data frame 
  par_temp_df <- par_sum_df_list[[sc]]
  
  #replace periods with underscores
  colnames(par_temp_df) <- gsub("\\.", "_", colnames(par_temp_df))
  
  #join maternal information with the parentage summary 
  par_temp_df <- left_join(par_temp_df, UHA_database, 
                           by=c('Mother_ID' = 'Tissue_ID'))
  
  #rename the species data frame to the maternal species 
  par_temp_df <- par_temp_df %>% rename("Maternal_Species" = "Species",
                                        "Maternal_Longitude" = "Longitude",
                                        "Maternal_Latitude" = "Latitude",
                                        "Maternal_Accession" = "Accession_Number")
  
  
  ##reorg data frame 
  #Adding in Species Information for Maternal and Paternal trees
  keep_col_ID <- c("Offspring_ID","Mother_ID", "Candidate_father_ID", 
                   "Maternal_Species", "Maternal_Longitude",
                   "Maternal_Latitude", "Maternal_Accession")
  
  # Narrowing the data again after the join of data sets
  par_temp_df <- par_temp_df[keep_col_ID] 
  
  #add paternal information 
  par_temp_df <- left_join(par_temp_df, UHA_database, 
                              by=c('Candidate_father_ID' = 'Tissue_ID'))  
  
  #rename columns 
  par_temp_df <- par_temp_df %>% rename('Candidate_Father_Species' = 'Species',
                              "Candidate_Father_Longitude" = "Longitude",
                              "Candidate_Father_Latitude" = "Latitude",
                              "Candidate_Father_Accession" = "Accession_Number")
  
  #reduce by empty columns 
  keep_col_ID2 <- c("Offspring_ID","Mother_ID", "Candidate_father_ID",
                    "Maternal_Species", "Maternal_Longitude", "Maternal_Latitude",
                    "Maternal_Accession", "Candidate_Father_Species", 
                    "Candidate_Father_Longitude", "Candidate_Father_Latitude",
                    "Candidate_Father_Accession")
  
  #reduce data frame by populated columns
  par_temp_df <- par_temp_df[keep_col_ID2]
  
  ###do half sibling analysis 
  ##Mikaely Evans' code 
  #initialize columns in the data frame 
  par_temp_df$Half_Sibs <- NA  # Made three new columns for this analysis
  par_temp_df$M_Accession_Abrv <- NA  
  par_temp_df$F_Accession_Abrv <- NA
  
  #Abbreviating the accession numbers for the maternal and paternal trees was necessary to assign half sibling status because it is only necessary to look at the first 6 characters in the accession number to know if the trees came from the same lineage.
  par_temp_df$M_Accession_Abrv <- substr(par_temp_df$Maternal_Accession, 0, 6)
  
  #This line adds the Maternal accession abbreviation to the M_Accession_Abrv column
  par_temp_df$F_Accession_Abrv <- substr(par_temp_df$Candidate_Father_Accession, 0, 6)
  
  #This line adds the paternal accession abbreviation to the F_Accession_Abrv column
  par_temp_df <- par_temp_df %>%
    mutate("Half_Sibs" = case_when(M_Accession_Abrv == F_Accession_Abrv ~ TRUE,
                                   M_Accession_Abrv != F_Accession_Abrv ~ FALSE))
  
  ##Add geographic information 
  #create a column for distance between mom and dad 
  par_temp_df$dist_par <- NA
  
  #loop to calculate distance between parents
  for(dist in 1:nrow(par_temp_df)){
    
    par_temp_df$dist_par[dist] <- distm(par_temp_df[dist,5:6], 
                                        par_temp_df[dist,9:10],
                                        fun=distGeo)
    
  }
  
  ###Code to add hybrid status 
  par_temp_df$Hybrid_Status <- NA
  par_temp_df <- par_temp_df %>%
    mutate(Hybrid_Status = case_when(Maternal_Species == Candidate_Father_Species ~ FALSE,
                                     Maternal_Species != Candidate_Father_Species ~ TRUE))
  
  #write to file 
  write.csv(par_temp_df, paste0("Data_Files/CSV_Files/UHA_",full_scen[[sc]], "_analysis_df.csv"),
            row.names = FALSE)
  
  
}


# ###################################
# #     Analyze Post Parentage      #
# ###################################
# #sum df 
# null_all_comp_df <- matrix(nrow = length(all_loc_par_sum$Candidate_father_ID),
#                            ncol = 3)
# #compare the two columns
# null_all_comp_df[,1] <- all_loc_par_sum$Candidate_father_ID == red_loc_par_sum$Candidate_father_ID #true is 1, false is 0
# 
# #add a column for all loci pair LOD score
# null_all_comp_df[,2] <- all_loc_par_sum$Pair_LOD_score
# 
# #add a column for red loci pair LOD score
# null_all_comp_df[,3] <- red_loc_par_sum$Pair_LOD_score
# 
# colnames(null_all_comp_df) <- c("Assigned_Father_Same", "All_Loc_LOD", "Red_Loc_LOD")
# rownames(null_all_comp_df) <- all_loc_par_sum$Offspring_ID
# 
# 
# 
# #subset by mismatch
# null_all_dif_df <- as.data.frame(null_all_comp_df[null_all_comp_df[,1] == FALSE,])
# 
# #add column greater
# null_all_dif_df$loc_greater <- NA
# 
# for(n in 1:length(null_all_dif_df[,1])){
#   if(null_all_dif_df[n,2] > null_all_dif_df[n,3]){
#     
#     null_all_dif_df$loc_greater[[n]] <- colnames(null_all_dif_df)[[2]]  
#     
#   }else{
#     null_all_dif_df$loc_greater[[n]] <- colnames(null_all_dif_df)[[3]]  
#   }
# }

# ############ NOT IN USE -----------
# #save as a data frame 
# #null_all_comp_df <- as.data.frame(null_all_comp_df)
# 
# #summarize - how many rows are false?
# #mismatch_names <- rownames(null_all_comp_df[null_all_comp_df[,1] == 0,])
# #mismatch_num <- length(null_all_comp_df[null_all_comp_df[,1] == 0,][,1])
# #15 individuals with 
# # #calculate mean distance between parents 
# # UHA_dist_matrix <- matrix(nrow = length(unique(full_parentage$Mother_ID)),
# #                           ncol = 1)
# # 
# # for(m in 1:length(unique(full_parentage$Mother_ID))){
# #   
# #   UHA_dist_matrix[m,1] <- mean(full_parentage[full_parentage$Mother_ID == unique(full_parentage$Mother_ID)[[m]],][,11])
# #   
# # }
