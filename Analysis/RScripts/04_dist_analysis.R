#####################
#     Libraries     #
#####################

library(tidyverse)
library(geosphere)

######################
#     Load files     #
######################

#set working directory
setwd("../..")

#load in the tissue database, remove offspring which have no coordinates
UHA_db <- read.csv("Data_Files/CSV_Files/ARCHIVED_USBG_Hybrid_Acorn_Tissue_Database.csv")

#read in cleaned data file
UHA_score_clean_df <- read.csv("Data_Files/CSV_Files/UHA_score_clean_df.csv")  

#load in parentage results 
parentage_results <- read.csv("Data_Files/CSV_Files/UHA_full_parentage.csv")

###################################
#     Reorganize data frames      #
###################################
#remove individuals with too much missing data
UHA_clean_df <- UHA_db[UHA_db$Tissue_ID %in% UHA_score_clean_df$Tissue_ID,]

#remove offspring from the data frame
UHA_par_df <- UHA_clean_df %>% 
                filter(!is.na(Longitude) | !is.na(Latitude))

#update parentage analysis df
UHA_clean_par_df <- parentage_results %>% 
                      rename(Father_ID = Candidate_father_ID) %>%
                            filter(!is.na(Father_ID)) 

#########################################
#     Prepare for Distance Analysis     #
#########################################

###Create data frames with distances
#Create a column of the crosses of every parent individual
all_potential_combo <- crossing(UHA_clean_df$Tissue_ID, 
                                UHA_clean_df$Tissue_ID)
colnames(all_potential_combo) <- c("Parent_1", "Parent_2")

#remove rows where the parents are the same individual (selfing)
potential_combo_dedup <- filter(all_potential_combo, Parent_1 != Parent_2)

#make columns to store distance between parents and parent species
potential_combo_dedup$dist <- NA 
potential_combo_dedup$Parent_1_species <- NA
potential_combo_dedup$Parent_2_species <- NA

#loop to calculate distance between parents 
for(d in 1:nrow(potential_combo_dedup)){
  Parent_1 <- potential_combo_dedup$Parent_1[d]
  Parent_2 <- potential_combo_dedup$Parent_2[d]
  
  #access the original tissue database via the parents
  Parent_1_Database_row <- filter(UHA_clean_df, Tissue_ID == Parent_1)
  Parent_2_Database_row <- filter(UHA_clean_df, Tissue_ID == Parent_2)
  
  #use distGeo to get distance in m between the lat long points of the 2 parents
  potential_combo_dedup$dist[d] <- distGeo(c(Parent_1_Database_row$Longitude, Parent_1_Database_row$Latitude),
                                           c(Parent_2_Database_row$Longitude, Parent_2_Database_row$Latitude)) 
  
  #record the species of both parents
  potential_combo_dedup$Parent_1_species[d] <- Parent_1_Database_row$Species 
  potential_combo_dedup$Parent_2_species[d] <- Parent_2_Database_row$Species
  
}



#reorganize the combination data frame to fix species naming conventions and to
#have a hybrid column 
potential_combo_info <- potential_combo_dedup %>%
  mutate(Parent_1_species = str_replace_all(Parent_1_species, "Quercus", "Q."), 
         Parent_2_species = str_replace_all(Parent_2_species, "Quercus", "Q.")) %>% #replace all instances of Quercus with Q. to ensure all species names are formatted the same
  mutate(Parental_species_match = case_when(Parent_1_species == Parent_2_species ~ "Conspecific",
                                            Parent_1_species != Parent_2_species ~ "Heterospecific")) #if the species of the two parents matches then they are conspecific, if not they are heterospecific

#list of maternal names
mom_IDs <- unique(parentage_results$Mother_ID) 
#filter combintation df by the possible mothers 
relevant_potential_combos <- potential_combo_info %>%
                                filter(Parent_1 %in% mom_IDs) 

#replaces relevant_parentage_results
#combine all columns for distance analysis 
par_results_df <- left_join(UHA_clean_par_df, 
                                        select(relevant_potential_combos, 
                                               c(Parent_1, Parent_2, dist)), 
                                                  join_by(Mother_ID == Parent_1,
                                                      Father_ID == Parent_2))

############### Create data frames
# right now just taking the mean dist of successful dads vs possible dads
#create df that has the mean distance of the 5 closest real fathers to each
#maternal tree (without ties, only possible with slice_min, using
#that instead of top_n)
rf_mean_small_df <- par_results_df %>%
                          group_by(Mother_ID) %>%
                            slice_min(dist, n =5, with_ties = FALSE) %>%
                              summarise(Mean_smallest_real_dists = mean(dist, 
                                                                        na.rm=TRUE))
#same as above but doing 5 farthest dists
rf_mean_large_df <- par_results_df %>% 
                      group_by(Mother_ID) %>%
                       slice_max(dist, n =5, with_ties = FALSE) %>%
                        summarise(Mean_largest_real_dists = mean(dist, na.rm=TRUE))


#create df that has the mean distance of the 5 closest possible 
#fathers to each maternal tree (including those with less than 5 offspring) 
#Note that I don't need to call distinct here because there is only 
#one entry per mother/father combo
pf_mean_small_df <- relevant_potential_combos %>%
                      group_by(Parent_1, Parental_species_match) %>%
                        slice_min(dist, n =5, with_ties = FALSE) %>%
                        summarise(Mean_smallest_potential_dists = mean(dist,
                                                                       na.rm=TRUE))

#same as abovebut with 5 farthest dists
pf_mean_large_df <- relevant_potential_combos %>%
                      group_by(Parent_1, Parental_species_match) %>%
                        slice_max(dist, n =5, with_ties = FALSE) %>%
                          summarise(Mean_largest_potential_dists = mean(dist, 
                                                                        na.rm=TRUE))

#################################
#     Run Distance Analysis     #
#################################

#summarize data by mean and min distance to real fathers and by proportion of 
#offspring that are hybrids (with heterospecific fathers)
rf_df <- par_results_df %>%
          group_by(Mother_ID) %>%
          summarise(Mean_real_dist = mean(dist, na.rm=TRUE), 
            Min_real_dist = min(dist, na.rm=TRUE), 
            Max_real_dist = max(dist, na.rm=TRUE), 
            Prop_hybrids = mean(Hybrid, na.rm = TRUE)) %>%
            left_join(., rf_mean_small_df, join_by(Mother_ID == Mother_ID)) %>% #add the Mean_smallest_dists data 
            left_join(., rf_mean_large_df, join_by(Mother_ID == Mother_ID)) #add the Mean_largest_dists data 

# Make dataset ("potential_fathers_summary") containing all of the desired 
#summarized categories from the possible combinations of mothers and fathers 
#data (mean distance of potential fathers, minimum distance of potential fathers, 
#mean distance of closest 5 potential fathers). 
#NOTE: this is NOT based off of the offspring data meaning that the 5 
#closest fathers will be unique and the mean distance to fathers will be 
#unweighted (because we have no way to know potential fathers relative success 
#(i.e. how many offspring they may produce)) summarize data from data 
#with all combos of possible conspecific and heterospefic individuals with 
#each real mom
pf_df <- relevant_potential_combos %>% 
           group_by(Parent_1, Parental_species_match) %>%
           summarise(Mean_potential_dist = mean(dist, na.rm=TRUE), 
            Min_potential_dist = min(dist, na.rm=TRUE), 
            Max_potential_dist = max(dist, na.rm=TRUE))  %>%
  left_join(., pf_mean_small_df, join_by(Parent_1 == Parent_1, Parental_species_match == Parental_species_match)) %>% #add the Mean_smallest_dists data
  left_join(., pf_mean_large_df, join_by(Parent_1 == Parent_1, Parental_species_match == Parental_species_match)) %>% #add the Mean_largest_dists data
  left_join(., select(rf_df, c(Mother_ID, Prop_hybrids)), join_by(Parent_1 == Mother_ID)) #add the proportion of hybrids had by each mother (from real_fathers_summary)

#now summarize this 
# The linear model for the mean distance of real fathers to a 
#given mother is significant with a positive slope indicating 
#that the higher the proportion of hybrid offspring a mother has, 
#the greater the mean realized pollination dsitance for that mother

###create a data frame to save the adjusted R-squared and p-value
hybrid_pop_df <- matrix(nrow = 1, ncol = 1)

#adjusted R-squared and 
hybrid_pop_df <- t(as.data.frame(c(summary(lm(formula = Prop_hybrids~Mean_real_dist, data=rf_df))[[9]],
                                   summary(lm(formula = Prop_hybrids~Mean_real_dist, data=rf_df))$coefficients[8])))

#add columns and row labels
rownames(hybrid_pop_df) <- "Case"
colnames(hybrid_pop_df) <- c("R2", "p-value")

#write out plot comparing real fathers plot
rf_df %>%
  ggplot() +
  geom_point(aes(y = Prop_hybrids, x = Mean_real_dist)) +
  geom_smooth(aes(y = Prop_hybrids, x = Mean_real_dist), method='lm', formula= y~x) +
  xlab("Mean Distance Between Parents (m)") + ylab("Proportion of Hybrid Offspring") +
  theme_classic()
