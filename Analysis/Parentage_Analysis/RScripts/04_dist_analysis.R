#####################
#     Libraries     #
#####################

library(tidyverse)
library(geosphere)

######################
#     Load files     #
######################

#set working directory
setwd("../../..")

#load in the tissue database, remove offspring which have no coordinates
UHA_database <- read.csv("Data_Files/CSV_Files/UHA_database_clean.csv")

#limit to one case right now 
genotype_df <- read.csv("Data_Files/CSV_Files/UHA_AL_par_sum.csv")

par_df <- read.csv("Results/Parentage_Results/CSV_Files/UHA_all_loci_analysis_df.csv")

###################################
#     Reorganize data frames      #
###################################

# #Loop over all data files 
# for(par in seq_along(par_scen_df)){
#   
#   temp_score_df <- score_df_list[[1]]
  
#remove offspring from the data frame - only modeling parents
UHA_par_df <- UHA_database %>% 
                filter(Parent_Offspring == "P")
  
#create data frame for offspring 
offspring_df <- par_df %>%
                    filter(!is.na(Candidate_father_ID))

###Create data frames with distances
#Create a column of the crosses of every parent individual
all_potential_combo <- crossing(UHA_par_df$Tissue_ID, 
                                  UHA_par_df$Tissue_ID)
#rename columns
colnames(all_potential_combo) <- c("Parent_1", "Parent_2")

#remove rows where the parents are the same individual (selfing)
potential_combo_dedup <- filter(all_potential_combo, Parent_1 != Parent_2)

#make columns to store distance between parents and parent species
potential_combo_dedup$dist <- NA 
potential_combo_dedup$Parent_1_species <- NA
potential_combo_dedup$Parent_2_species <- NA

#loop to calculate distance between parents 
for(d in 1:nrow(potential_combo_dedup)){
  
  #initialize first and second parent
  Parent_1 <- potential_combo_dedup$Parent_1[d]
  Parent_2 <- potential_combo_dedup$Parent_2[d]
  
  #access the original tissue database via the parents
  Parent_1_Database_row <- filter(UHA_database, Tissue_ID == Parent_1)
  Parent_2_Database_row <- filter(UHA_database, Tissue_ID == Parent_2)
  
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
  mutate(Parent_1_species = str_replace_all(Parent_1_species, "Q.", "Quercus"), 
         Parent_2_species = str_replace_all(Parent_2_species, "Q.", "Quercus")) %>% #replace all instances of Q.
  mutate(Parental_species_match = case_when(Parent_1_species == Parent_2_species ~ "Conspecific",
                                            Parent_1_species != Parent_2_species ~ "Heterospecific")) #if the species of the two parents matches then they are conspecific, if not they are heterospecific

#list of maternal names
mom_IDs <- unique(par_df$Maternal_ID) 
#filter combintation df by the possible mothers 
relevant_potential_combos <- potential_combo_info %>%
                                filter(Parent_1 %in% mom_IDs) 

#replaces relevant_parentage_results
#combine all columns for distance analysis 
par_results_df <- left_join(par_df, select(relevant_potential_combos, 
                                               c(Parent_1, Parent_2, dist)), 
                                                  join_by(Maternal_ID == Parent_1,
                                                          Candidate_father_ID == Parent_2))

############### Create data frames
# right now just taking the mean dist of successful dads vs possible dads
#create df that has the mean distance of the 5 closest real fathers to each
#maternal tree (without ties, only possible with slice_min, using
#that instead of top_n)
rf_mean_small_df <- par_results_df %>%
                          group_by(Maternal_ID) %>%
                            slice_min(dist, n =5, with_ties = FALSE) %>%
                              summarise(Mean_smallest_real_dists = mean(dist, 
                                                                        na.rm=TRUE))
#same as above but doing 5 farthest dists
rf_mean_large_df <- par_results_df %>% 
                      group_by(Maternal_ID) %>%
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
          group_by(Maternal_ID) %>%
          summarise(Mean_real_dist = mean(dist, na.rm=TRUE), 
            Min_real_dist = min(dist, na.rm=TRUE), 
            Max_real_dist = max(dist, na.rm=TRUE), 
            Prop_hybrids = mean(Hybrid_Status, na.rm = TRUE)) %>%
            left_join(., rf_mean_small_df, join_by(Maternal_ID == Maternal_ID)) %>% #add the Mean_smallest_dists data 
            left_join(., rf_mean_large_df, join_by(Maternal_ID == Maternal_ID)) #add the Mean_largest_dists data 

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
  left_join(., select(rf_df, c(Maternal_ID, Prop_hybrids)), join_by(Parent_1 == Maternal_ID)) #add the proportion of hybrids had by each mother (from real_fathers_summary)

#now summarize this 
# The linear model for the mean distance of real fathers to a 
#given mother is significant with a positive slope indicating 
#that the higher the proportion of hybrid offspring a mother has, 
#the greater the mean realized pollination dsitance for that mother

###create a data frame to save the adjusted R-squared and p-value
rf_hybrid_prop_dist_df <- matrix(nrow = 5, ncol = 2)

# #adjusted R-squared and 
# hybrid_pop_df <- t(as.data.frame(c(summary(lm(formula = Prop_hybrids~Mean_real_dist, data=rf_df))[[9]],
#                                    summary(lm(formula = Prop_hybrids~Mean_real_dist, data=rf_df))$coefficients[8])))

#save p-values and r-squared results
rf_hybrid_prop_dist_df[1,1] <- summary(lm(formula = Prop_hybrids~Mean_real_dist, data=rf_df))$adj.r.squared
rf_hybrid_prop_dist_df[1,2] <- summary(lm(formula = Prop_hybrids~Mean_real_dist, data=rf_df))$coefficients[2,4]
rf_hybrid_prop_dist_df[2,1] <- summary(lm(formula = Prop_hybrids~Min_real_dist, data=rf_df))$adj.r.squared
rf_hybrid_prop_dist_df[2,2] <- summary(lm(formula = Prop_hybrids~Min_real_dist, data=rf_df))$coefficients[2,4]
rf_hybrid_prop_dist_df[3,1] <- summary(lm(formula = Prop_hybrids~Max_real_dist, data=rf_df))$adj.r.squared
rf_hybrid_prop_dist_df[3,2] <- summary(lm(formula = Prop_hybrids~Max_real_dist, data=rf_df))$coefficients[2,4]
rf_hybrid_prop_dist_df[4,1] <- summary(lm(formula = Prop_hybrids~Mean_smallest_real_dists, data=rf_df))$adj.r.squared
rf_hybrid_prop_dist_df[4,2] <- summary(lm(formula = Prop_hybrids~Mean_smallest_real_dists, data=rf_df))$coefficients[2,4]
rf_hybrid_prop_dist_df[5,1] <- summary(lm(formula = Prop_hybrids~Mean_largest_real_dists, data=rf_df))$adj.r.squared
rf_hybrid_prop_dist_df[5,2] <- summary(lm(formula = Prop_hybrids~Mean_largest_real_dists, data=rf_df))$coefficients[2,4]

#add columns and row labels
rownames(rf_hybrid_prop_dist_df) <- c("Mean_Par_Dist", "Min_Par_Dist", "Max_Par_Dist",
                              "Mean_Smallest_Dist", "Mean_Largest_Dist")
colnames(rf_hybrid_prop_dist_df) <- c("R2", "p-value")

#write out plot comparing real fathers plot
png("Results/Pairwise_Distance_Analysis/rf_dist_par_prop.png", width = 900,
    height = 600)
rf_df %>%
  ggplot() +
  geom_point(aes(y = Prop_hybrids, x = Mean_real_dist)) +
  geom_smooth(aes(y = Prop_hybrids, x = Mean_real_dist), method='lm', formula= y~x) +
  xlab("Mean Distance Between Parents (m)") + ylab("Proportion of Hybrid Offspring") +
  theme_classic()
dev.off()

# The linear model for the mean distance of the shortest 5 distances of real fathers
#to a given mother shows no relationship between the distance of the closest 
#5 realized fathers to a given mother and the proportion of offspring that were 
#hybrids for that mother.
rf_df %>%
  ggplot() +
  geom_point(aes(y = Prop_hybrids, x = Mean_smallest_real_dists)) +
  geom_smooth(aes(y = Prop_hybrids, x = Mean_smallest_real_dists), 
              method='lm', formula= y~x) +
  xlab("Mean Distance Between Parents (m) - Shortest Distances") + 
  ylab("Proportion of Hybrid Offspring") +
  theme_classic()

# The linear model for the mean distance of the largest 5 distances of real fathers
#to a given mother 
rf_df %>%
  #filter(Mean_largest_real_dists < 500) %>% #to remove the largest value individual if you want
  ggplot() +
  geom_point(aes(y = Prop_hybrids, x = Mean_largest_real_dists)) +
  geom_smooth(aes(y = Prop_hybrids, x = Mean_largest_real_dists), 
              method='lm', formula= y~x) +
  xlab("Mean Distance Between Parents (m) - Largest Distances") + 
  ylab("Proportion of Hybrid Offspring") +
  theme_classic()

# note that when you remove the mom with the largest values (> 500), 
#still significant
summary(lm(formula = Prop_hybrids~Mean_largest_real_dists, 
           data=filter(rf_df, Mean_largest_real_dists < 500)))

write.csv(rf_hybrid_prop_dist_df, "Results/Pairwise_Distance_Analysis/rf_hybrid_prop_dist.csv")

#####Potential fathers to a given mother df 
###create a data frame to save the adjusted R-squared and p-value
pf_hybrid_prop_dist <- matrix(nrow = 6, ncol = 2)

# The linear model for the mean distance of a potential fathers to a given mother split by conspecific and heterospecific fathers. There is no relationship between the proportion of hybrid offspring and the mean distances of all potential het or con fathers.
pf_hybrid_prop_dist[1,1] <- summary(lm(formula = Prop_hybrids~Mean_potential_dist, 
                                       data=filter(pf_df, 
                                                   Parental_species_match == "Conspecific")))$adj.r.squared
pf_hybrid_prop_dist[1,2] <- summary(lm(formula = Prop_hybrids~Mean_potential_dist, 
                                       data=filter(pf_df, 
                                                   Parental_species_match == "Conspecific")))$coefficients[2,4]
#mean distance to potential fathers, conspecific 
pf_hybrid_prop_dist[2,1] <- summary(lm(formula = Prop_hybrids~Mean_potential_dist, data=filter(pf_df, Parental_species_match == "Heterospecific")))$adj.r.squared
pf_hybrid_prop_dist[2,2] <- summary(lm(formula = Prop_hybrids~Mean_potential_dist, data=filter(pf_df, Parental_species_match == "Heterospecific")))$coefficients[2,4]


#The linear model for the smallest distance of a potential father to a given 
#mother split by conspecific and heterospecific fathers. 
#There is no relationship between the proportion of hybrid offspring and the 
#minimum distances of all potential het or con fathers.
pf_hybrid_prop_dist[3,1] <- summary(lm(formula = Prop_hybrids~Min_potential_dist, data=filter(pf_df, Parental_species_match == "Conspecific")))$adj.r.squared
pf_hybrid_prop_dist[3,2]<- summary(lm(formula = Prop_hybrids~Min_potential_dist, data=filter(pf_df, Parental_species_match == "Conspecific")))$coefficients[2,4]


pf_hybrid_prop_dist[4,1] <- summary(lm(formula = Prop_hybrids~Min_potential_dist, data=filter(pf_df, Parental_species_match == "Heterospecific")))$adj.r.squared
pf_hybrid_prop_dist[4,2]<- summary(lm(formula = Prop_hybrids~Min_potential_dist, data=filter(pf_df, Parental_species_match == "Heterospecific")))$coefficients[2,4]


# The linear model for the mean distance of the shortest 5 distances of potential fathers to a given mother split by conspecific and heterospecific fathers. There is no relationship between the proportion of hybrid offspring and the 5 smallest distances of all potential het or con fathers.
pf_hybrid_prop_dist[5,1] <- summary(lm(formula = Prop_hybrids~Mean_smallest_potential_dists, data=filter(pf_df, Parental_species_match == "Conspecific")))$adj.r.squared
pf_hybrid_prop_dist[5,2] <- summary(lm(formula = Prop_hybrids~Mean_smallest_potential_dists, data=filter(pf_df, Parental_species_match == "Conspecific")))$coefficients[2,4]

pf_hybrid_prop_dist[6,1] <- summary(lm(formula = Prop_hybrids~Mean_smallest_potential_dists, data=filter(pf_df, Parental_species_match == "Heterospecific")))$adj.r.squared
pf_hybrid_prop_dist[6,2] <- summary(lm(formula = Prop_hybrids~Mean_smallest_potential_dists, data=filter(pf_df, Parental_species_match == "Heterospecific")))$coefficients[2,4]

#add rownames and colnames
rownames(pf_hybrid_prop_dist) <- c("Mean_Dist_Con", "Mean_Dist_Het", "Min_Dist_Con",
                                   "Min_Dist_Het", "Max_Dist_Con", "Max_Dist_Het")
colnames(pf_hybrid_prop_dist) <- c("R2", "p-value")

write.csv(pf_hybrid_prop_dist, "Results/Pairwise_Distance_Analysis/pf_hybrid_prop_dist.csv")

#########################################
#     Resampling Distance Analysis      #
#########################################

# A for loop that, across each mom, randomly samples the pool of the 
#possible conspecific dads (of which there are 20 out of the 280 adult 
#individuals we have data on) 1000 times 
#and records the distance between that dad and the given mom 


set.seed(2024) #since we are using the sample function, set the seed for repeatability 

#create the table that will hold the results of the for loop with 3 columns: Mom, Dad, and dist
parents_dists_table_full <- tibble(Mom = character(), Dad = character(), dist = character())


for(i in 1:length(unique(relevant_potential_combos$Parent_1))){
  parents_dists_table_small <- tibble(Mom = character(), Dad = character(), dist = character()) #make a small version of the final table that will be overwritten with each new mom
  mom <- unique(relevant_potential_combos$Parent_1)[i]  #get the ID of the given mom
  mom_specific_combos <- relevant_potential_combos %>%
    filter(Parent_1 == mom) %>%  #filter all possible combos to only those with the given mom as Parent_1
    filter(Parental_species_match == "Conspecific")  #filter all possible combos with the given mom to only those that are the same species as the given mom
  
  possible_dads <- mom_specific_combos$Parent_2 #make a vector with a list of all of the possible dads by pulling the Parent_2 column of the mom_specific_combos df
  
  dads <- sample(possible_dads, size=1000, replace=T)  #randomly draw 1000 samples from all possible dads with replacement
  
  #a for loop that will loop through each of the 1000 sampled dads will make and add a row for the parents_dists_table_small 
  for(x in 1:length(dads)){
    dad <- dads[x] #get the ID of the given dad
    parent_combo_row_for_table <- tibble(Mom = NA, Dad = NA, dist = NA) #creating table that will be a single row and will be overwritten during each iteration of the loop
    parent_combo_row_for_table$Mom <- mom #add ID of the given mom to table
    parent_combo_row_for_table$Dad <- dad #add ID of the given dad to table
    parent_combo_row_for_table$dist <- filter(mom_specific_combos, Parent_1 == mom & Parent_2 == dad)$dist #add dist from mom to dad in table (by matching the mom and dad's IDs to the mom_specific_combos df)
    
    parents_dists_table_small <- rbind(parents_dists_table_small, parent_combo_row_for_table) #bind the row of the given dad to the  table of the given mom
  }
  parents_dists_table_full <- rbind(parents_dists_table_full, parents_dists_table_small)  #bind the table of the given given to the  table of all the moms
}


# Make a plot that compares the real distances between parents 
#(from relevant_parentage_results) to the potential distance between
#parents (from parents_dists_table_full)

parentage_results_for_hists <- par_results_df %>%
  rename(Mom = Maternal_ID, 
         Dad = Candidate_father_ID) %>%
  select(Mom, Dad, dist, Hybrid_Status) %>%
  filter(Hybrid_Status == F) %>%  #we are only looking at distance between conspecific pairs in this analysis (no hydrbid offspring)
  filter(!is.na(dist))  #remove entries where dist is NA because some trees that we have tissue from do not have long/lat values

#Create the a histogram that shows the difference in the real and theoretical distances between parents 

dist_hist <- ggplot() +
  geom_histogram(data = parents_dists_table_full, 
                 aes(x = dist, y = ..density.., fill = "theoretical"), 
                 alpha = .5) +
  geom_histogram(data = parentage_results_for_hists, 
                 aes(x = dist, y = ..density.., fill = "real"), 
                 alpha = .5) +
  labs(fill = "data set") + 
  scale_fill_manual(values = c("theoretical" = "red", "real" = "blue")) + 
  #facet_wrap(~Mom) +  # facet_wrap by Mom to see how each individual Maternal tree did compared to bootstrap
  theme_classic()

##write out histogram 
png("Results/Pairwise_Distance_Analysis/sim_real_dist_hist.png", width = 1000, height = 850)

dist_hist

dev.off()

#compare relationships 
ks.test(x = parentage_results_for_hists$dist, 
        y = parents_dists_table_full$dist)

