###This script creates barplots of the plots presented in the manuscript
##This code was developed and tested by Mikaely Evans in 2023 


#####################
#     Libraries     #
#####################

library(tidyverse)
library(ggplot2)

###########################
#     Load Data Files     #
###########################

#set working directory
setwd("../../..")

#load parentage result dfs
par_scen_df_list <- list.files(path = "Data_Files/CSV_Files/",pattern = "analysis_df.csv")

#reorder 
par_scen_df_list <- list(par_scen_df_list[[1]], par_scen_df_list[[2]],
                         par_scen_df_list[[4]], par_scen_df_list[[3]])

#load full scenario data frames 
full_scen <- c("all_loci_AF", #all loci included with all father assignments 
               "all_loci_HCF", #all loci with only high confidence fathers included
               "red_loci_AF", #reduced loci with all father assignments
               "red_loci_HCF" #reduce loci with only high confidence father assignments included
               )

#read in these files as CSVs
par_scen_df <- list()

for(df in 1:length(par_scen_df_list)){
  
  par_scen_df[[df]] <- read.csv(paste0("Data_Files/CSV_Files/", par_scen_df_list[[df]]))
  
}

###########################
#     Visualizations      #
###########################

####barplots of candidate fathers for each mother, over scenario 

#plot bar graphs of the offspring count for each mother and which individual is the father
png(paste0("Results/Parentage_Figures/", full_scen[[1]], "_mat_off.png"), 
    width = 5000, height = 3500, res = 600)
par_scen_df[[1]] %>%
  ggplot() +
  geom_bar(aes(y = sort(Candidate_father_ID))) +
  facet_wrap(~`Mother_ID`) + 
  scale_x_continuous(n.breaks = 9) +
  labs(title = "Count of Offspring per Candidate Father and Maternal Tree Pairs", 
       y = "Candidate Father ID", x = "Count of Offspring")
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[2]], "_mat_off.png"), 
    width = 5000, height = 3500, res = 600)
par_scen_df[[2]] %>%
  ggplot() +
  geom_bar(aes(y = sort(Candidate_father_ID))) +
  facet_wrap(~`Mother_ID`) + 
  scale_x_continuous(n.breaks = 9) +
  labs(title = "Count of Offspring per Candidate Father and Maternal Tree Pairs", 
       y = "Candidate Father ID", x = "Count of Offspring")
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[3]], "_mat_off.png"), 
    width = 5000, height = 3500, res = 600)
par_scen_df[[3]] %>%
  ggplot() +
  geom_bar(aes(y = sort(Candidate_father_ID))) +
  facet_wrap(~`Mother_ID`) + 
  scale_x_continuous(n.breaks = 9) +
  labs(title = "Count of Offspring per Candidate Father and Maternal Tree Pairs", 
       y = "Candidate Father ID", x = "Count of Offspring")
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[4]], "_mat_off.png"), 
    width = 5000, height = 3500, res = 600)
par_scen_df[[4]] %>%
  ggplot() +
  geom_bar(aes(y = sort(Candidate_father_ID))) +
  facet_wrap(~`Mother_ID`) + 
  scale_x_continuous(n.breaks = 9) +
  labs(title = "Count of Offspring per Candidate Father and Maternal Tree Pairs", 
       y = "Candidate Father ID", x = "Count of Offspring")
dev.off()

###### Boxplot of distances between parents 
png(paste0("Results/Parentage_Figures/", full_scen[[1]], "_dist_par.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[1]] %>%
  group_by(c(Mother_ID)) %>% # 
  ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), y = dist_par)) +  
  expand_limits(y = c(0, 650)) +  # set limits for graph
  #theme_minimal() +  # set theme
  theme_bw() +  # set theme
  geom_boxplot(fill="darkolivegreen4") +
  geom_jitter(color = "grey", fill = "black", width = 0.3) +
  geom_text(data = . %>% count(Mother_ID), aes(label = n, y = 645), vjust = -0.5) + 
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)")

dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[2]], "_dist_par.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[2]] %>%
  group_by(c(Mother_ID)) %>% # 
  ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), y = dist_par)) +  
  expand_limits(y = c(0, 650)) +  # set limits for graph
  #theme_minimal() +  # set theme
  theme_bw() +  # set theme
  geom_boxplot(fill="darkolivegreen4") +
  geom_jitter(color = "grey", fill = "black", width = 0.3) +
  geom_text(data = . %>% count(Mother_ID), aes(label = n, y = 645), vjust = -0.5) + 
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)")

dev.off()


png(paste0("Results/Parentage_Figures/", full_scen[[3]], "_dist_par.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[3]] %>%
  group_by(c(Mother_ID)) %>% # 
  ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), y = dist_par)) +  
  expand_limits(y = c(0, 650)) +  # set limits for graph
  #theme_minimal() +  # set theme
  theme_bw() +  # set theme
  geom_boxplot(fill="darkolivegreen4") +
  geom_jitter(color = "grey", fill = "black", width = 0.3) +
  geom_text(data = . %>% count(Mother_ID), aes(label = n, y = 645), vjust = -0.5) + 
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)")

dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[4]], "_dist_par.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[4]] %>%
  group_by(c(Mother_ID)) %>% # 
  ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), y = dist_par)) +  
  expand_limits(y = c(0, 650)) +  # set limits for graph
  #theme_minimal() +  # set theme
  theme_bw() +  # set theme
  geom_boxplot(fill="darkolivegreen4") +
  geom_jitter(color = "grey", fill = "black", width = 0.3) +
  geom_text(data = . %>% count(Mother_ID), aes(label = n, y = 645), vjust = -0.5) + 
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)")

dev.off()

#loop over for all cases 

for(scen in 1:length(full_scen)){
  par_scen_df[[scen]] <- par_scen_df[[scen]] %>%
  mutate(`Parents Are` = case_when(Half_Sibs == FALSE ~ "Not Half Siblings",
                                                TRUE ~ "Half Siblings", 
                                                NA ~ "No Accession Number")
  )
}

#graph of half-sibling matings group by maternal ID
png(paste0("Results/Parentage_Figures/", full_scen[[1]], "_halfsib_dist.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[1]] %>%
  group_by(`Parents Are`) %>%  # group by half siblings to compare the status
  ggplot(aes(x = Mother_ID, y = dist_par, color = `Parents Are`)) +  
  geom_jitter(width = 0.2) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  scale_color_manual(values = c("cadetblue", "navy")) +
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw()
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[2]], "_halfsib_dist.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[2]] %>%
  group_by(`Parents Are`) %>%  # group by half siblings to compare the status
  ggplot(aes(x = Mother_ID, y = dist_par, color = `Parents Are`)) +  
  geom_jitter(width = 0.2) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  scale_color_manual(values = c("cadetblue", "navy")) +
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw()
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[3]], "_halfsib_dist.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[3]] %>%
  group_by(`Parents Are`) %>%  # group by half siblings to compare the status
  ggplot(aes(x = Mother_ID, y = dist_par, color = `Parents Are`)) +  
  geom_jitter(width = 0.2) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  scale_color_manual(values = c("cadetblue", "navy")) +
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw()
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[4]], "_halfsib_dist.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[4]] %>%
  group_by(`Parents Are`) %>%  # group by half siblings to compare the status
  ggplot(aes(x = Mother_ID, y = dist_par, color = `Parents Are`)) +  
  geom_jitter(width = 0.2) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  scale_color_manual(values = c("cadetblue", "navy")) +
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw()
dev.off()

#loop to add hybrid status to data frame 
for(scen in 1:length(full_scen)){
  
  par_scen_df[[scen]] <- par_scen_df[[scen]] %>%
    mutate(`Offspring Hybrid Status` = case_when(Hybrid_Status == FALSE ~ "Not Hybrid",
                                                 TRUE ~ "Hybrid"))
  
}

#####barplots of the percentage of hybrids 
png(paste0("Results/Parentage_Figures/", full_scen[[1]],"_hybrid_boxplot.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[1]] %>%
ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), 
                                    y = dist_par, 
                                  fill = `Offspring Hybrid Status`)) +  
  geom_boxplot() +
  scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw() 

dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[2]],"_hybrid_boxplot.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[2]] %>%
  ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), 
             y = dist_par, 
             fill = `Offspring Hybrid Status`)) +  
  geom_boxplot() +
  scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw() 

dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[3]],"_hybrid_boxplot.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[3]] %>%
  ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), 
             y = dist_par, 
             fill = `Offspring Hybrid Status`)) +  
  geom_boxplot() +
  scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw() 

dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[4]],"_hybrid_boxplot.png"), 
    res = 600, width = 5000, height = 3500)
par_scen_df[[4]] %>%
  ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), 
             y = dist_par, 
             fill = `Offspring Hybrid Status`)) +  
  geom_boxplot() +
  scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw() 

dev.off()

##create a list for 
#table to present the candidate fathers 
species_count_list <- list()

#loop to create species count lists for parent assignments
for(scen in 1:length(full_scen)){
  
  #create a data frame 
  species_count_list[[scen]] <- as.data.frame(table(par_scen_df[[scen]]$Candidate_Father_Species))
  
  #rename columns 
  names(species_count_list[[scen]]) <- c("Species", "Count")
  
  #organize data frame 
  species_count_list[[scen]]$Species <- factor(species_count_list[[scen]]$Species, 
                                  levels=species_count_list[[scen]]$Species[order(-species_count_list[[scen]]$Count)])
  
  
}


#now plot species count table
png(paste0("Results/Parentage_Figures/", full_scen[[1]], "_hybrid_barplot.png"),
    res = 600, width = 5000, height = 3500)
species_count_list[[1]] %>%
  ggplot(aes(x = Species, y=Count))+
  geom_bar(stat = "identity", fill = "darkgreen") + 
  labs(title="Count of Candidate Father Trees per Species", 
       x="Candidate Father Species") +
  theme_bw()
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[2]], "_hybrid_barplot.png"),
    res = 600, width = 5000, height = 3500)
species_count_list[[2]] %>%
  ggplot(aes(x = Species, y=Count))+
  geom_bar(stat = "identity", fill = "darkgreen") + 
  labs(title="Count of Candidate Father Trees per Species", 
       x="Candidate Father Species") +
  theme_bw()
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[3]], "_hybrid_barplot.png"),
    res = 600, width = 5000, height = 3500)
species_count_list[[3]] %>%
  ggplot(aes(x = Species, y=Count))+
  geom_bar(stat = "identity", fill = "darkgreen") + 
  labs(title="Count of Candidate Father Trees per Species", 
       x="Candidate Father Species") +
  theme_bw()
dev.off()

png(paste0("Results/Parentage_Figures/", full_scen[[4]], "_hybrid_barplot.png"),
    res = 600, width = 5000, height = 3500)
species_count_list[[4]] %>%
  ggplot(aes(x = Species, y=Count))+
  geom_bar(stat = "identity", fill = "darkgreen") + 
  labs(title="Count of Candidate Father Trees per Species", 
       x="Candidate Father Species") +
  theme_bw()
dev.off()

################## Unused Loops ------------- 
for(scen in seq_along(full_scen)){
  
  mat_off <- par_scen_df[[scen]] %>%
    ggplot() +
    geom_bar(aes(y = sort(Candidate_father_ID))) +
    facet_wrap(~`Mother_ID`) + 
    scale_x_continuous(n.breaks = 9) +
    labs(title = "Count of Offspring per Candidate Father and Maternal Tree Pairs", 
         y = "Candidate Father ID", x = "Count of Offspring")
  
  png(paste0("Results/Parentage_Figures/", full_scen[[scen]], "_mat_off.png"), 
      width = 5000, height = 3500, res = 600)
  
  mat_off
  
}
dev.off()


