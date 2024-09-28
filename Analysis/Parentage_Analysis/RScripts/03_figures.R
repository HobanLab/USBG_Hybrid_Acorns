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

#read in summary df for final figures
UHA_res_df <- read.csv("Results/Parentage_Results/CSV_Files/UHA_HCF_all_loci_analysis_df.csv")

#load in list of parentage summary results CSVs
par_sum_df <- list.files(path = "Results/Parentage_Results/CSV_Files",
                         pattern = "par_sum.csv")

#reorg list to fit order 
par_sum_df <- list(par_sum_df[[2]], par_sum_df[[1]],
                   par_sum_df[[4]], par_sum_df[[3]])

#load full scenario data frames 
full_scen <- c("all_loci_AF", #all loci included with all father assignments 
               "all_loci_HCF", #all loci with only high confidence fathers included
               "red_loci_AF", #reduced loci with all father assignments
               "red_loci_HCF" #reduce loci with only high confidence father assignments included
)


###################################################
#     Sumamry of Non Exculsion Probabilities      #
###################################################

##first, create a data frame to compare non-exclusion probablities
#create matrix to store non-exclusion probabilities
exc_prob_df <- matrix(nrow = length(par_sum_df),
                      ncol = 2)

#loop over all parentage summary data frames 
for(par in seq_along(par_sum_df)){
  
  #load in par sum_df
  temp_df <- read.csv(paste0("Results/Parentage_Results/CSV_Files/", 
                                     par_sum_df[[par]]))
  #rename colnames
  colnames(temp_df) <- gsub("\\.","_", colnames(temp_df))
  
  #now sum the columns 
  exc_prob_df[par,1] <- mean(temp_df[["First_parent_non_exclusion_probability"]])
  exc_prob_df[par,2] <- mean(temp_df[["Second_parent_non_exclusion_probability"]])
  
  
  
}

#name columns and rows 
rownames(exc_prob_df) <- c("all_loci", "all_loci_HCF", "red_loci", "red_loci_HCF")
colnames(exc_prob_df) <- c("First_parent_non_exclusion_probability",
                           "Second_parent_non_exclusion_probability")

#write out table 
write.csv(exc_prob_df, "./Results/Parentage_Results/CSV_Files/non_exclusion_probabilities.csv")

########################################
#     Summary Stat Table Creation      #
########################################

#loop over all parentage scenarios with a summary data file 
for(df in seq_along(par_scen_df_list)){
  
  ## create a summary data frame with the maternal IDs
  par_sum_stat_df <- matrix(nrow = length(mat_ids),
                            ncol = 11)
  
  #add the maternal ids in a column
  par_sum_stat_df[,1] <- mat_ids
  
  #loop over each maternal individual to create data frame 
  for(mat in seq_along(mat_ids)){
    
    #add a column for the offspring counts by mom 
    par_sum_stat_df[mat,2] <- length(par_scen_df[[df]][par_scen_df[[df]]$MT_ID == mat_ids[[mat]],"MT_ID"])
    
    #add a column for number of dads per mom 
    par_sum_stat_df[mat,3] <- length(unique(par_scen_df[[df]][par_scen_df[[df]]$MT_ID == mat_ids[[mat]],]$Candidate_father_ID))
    
    #add columns with the number of hybrid per mom and the number of hybrid fathers
    par_sum_stat_df[mat,4] <- length(par_scen_df[[df]][(par_scen_df[[df]]$MT_ID == mat_ids[[mat]]) & (par_scen_df[[df]]$Hybrid_Status == TRUE),]$Candidate_Father_Species)
    par_sum_stat_df[mat,5] <- length(unique(par_scen_df[[df]][(par_scen_df[[df]]$MT_ID == mat_ids[[mat]]) & (par_scen_df[[df]]$Hybrid_Status == TRUE),]$Candidate_Father_Species))
  
    #add half-sib columns 
    par_sum_stat_df[mat,6] <- length(par_scen_df[[df]][(par_scen_df[[df]]$MT_ID == mat_ids[[mat]]) & (par_scen_df[[df]]$Half_Sibs == TRUE),]$Offspring_ID)
    par_sum_stat_df[mat,7] <- length(unique(par_scen_df[[df]][(par_scen_df[[df]]$MT_ID == mat_ids[[mat]]) & (par_scen_df[[df]]$Half_Sibs == TRUE),]$Candidate_father_ID))
    
      if(par_sum_stat_df[mat,2] == 0){
        par_sum_stat_df[mat,8] <- 0
        par_sum_stat_df[mat,9] <- 0
        par_sum_stat_df[mat,10] <- 0
        par_sum_stat_df[mat,11] <- 0
        
      }else{
      
        par_sum_stat_df[mat,8] <- min(par_scen_df[[df]][(par_scen_df[[df]]$MT_ID == mat_ids[[mat]]),]$dist_par)
        par_sum_stat_df[mat,9] <- max(par_scen_df[[df]][(par_scen_df[[df]]$MT_ID == mat_ids[[mat]]),]$dist_par)
        
        #recode hybrid status of min distance 
        par_sum_stat_df[mat,10] <- length(par_scen_df[[df]][par_scen_df[[df]]$dist_par == min(par_scen_df[[df]][(par_scen_df[[df]]$MT_ID == mat_ids[[mat]]),]$dist_par) & par_scen_df[[df]]$Hybrid_Status == TRUE,]$Candidate_Father_Species)
        #recode hybrid status of min distance 
        par_sum_stat_df[mat,11] <- length(par_scen_df[[df]][par_scen_df[[df]]$dist_par == max(par_scen_df[[df]][(par_scen_df[[df]]$MT_ID == mat_ids[[mat]]),]$dist_par) & par_scen_df[[df]]$Hybrid_Status == TRUE,]$Candidate_Father_Species)
        
    }
    #add min and max dist columns for parents
    
    
    
    # #add colnames
    colnames(par_sum_stat_df) <- c("MT_ID", "Off_N",
                                  "Fathers_N", "Hybrid_Off_N",
                                  "Hybrid_Father_N", "Half_Sib_Off_N",
                                  "Half_Sib_Father_N", "Min_Dist", "Max_Dist",
                                  "Hybrid_Min", "Hybrid_Max")
    # #write out summary data frame with scenario 
    write.csv(par_sum_stat_df, paste0("Results/Parentage_Results/", full_scen[[df]], '_par_sum_stat_df.csv'),
             row.names = FALSE)
  }
  
  
  
}

###################
#     Figures     #
###################

###### Boxplot of distances between parents ------------------
png(paste0("Results/Parentage_Results/Figures/AL_HCF_dist_par.png"),
    res = 600, width = 5200, height = 3500)
UHA_res_df %>%
  group_by(c(Maternal_ID)) %>% # 
  ggplot(aes(x = fct_rev(fct_infreq(Maternal_ID)), y = dist_par)) +  
  expand_limits(y = c(0, 675)) +  # set limits for graph
  #theme_minimal() +  # set theme
  theme_bw() +  # set theme
  geom_boxplot(fill="darkolivegreen4", outlier.shape = NA) + # set color and remove outliers
  geom_jitter(aes(fill = Hybrid_Status), width = 0.2, size = 3.25, shape = 21, color = "black") +
  geom_text(data = . %>% count(Maternal_ID), aes(label = paste("n =", n), y = 665), vjust = -0.5) + 
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") + 
  scale_fill_manual(values = c("TRUE" = "hotpink", "FALSE" = "grey"),
                    labels = c("TRUE" = "Hybrid", "FALSE" = "Not a hybrid")) + # set color and titles for Hybrid Status
  labs(fill = "Offspring is: ", title = "Fig. 2: Distribution of Mating Distances Between Maternal and Paternal Trees") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.3))  # center the title

dev.off()

###### Boxplot of distances and half siblings ------------------

#replace column text for halfsibs
UHA_res_df[["Parents Are"]] <- dplyr::case_when(
  UHA_res_df$Half_Sibs == TRUE ~ "Half Siblings",
  UHA_res_df$Half_Sibs == FALSE ~ "Not Half Siblings"
)

#graph of half-sibling matings group by maternal ID
png(paste0("Results/Parentage_Results/Figures/AL_HCF_dist_par_halfsib.png"),
    res = 600, width = 5200, height = 3500)
UHA_res_df %>%
  group_by(`Parents Are`) %>%  # group by half siblings to compare the status
  ggplot(aes(x = fct_rev(fct_infreq(Maternal_ID)), y = dist_par, fill = `Parents Are`)) +  
  geom_jitter(aes(fill = `Parents Are`), width = 0.2, size = 3, shape = 21, color = "black") +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  scale_fill_manual(values = c("cadetblue", "navy")) +
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  labs(title = "Fig. 3: Distribution of Mating Distances Between Half-Sibling Parents") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.3))
dev.off()

# ###### Boxplot of distances between parents 
# png(paste0("Results/Parentage_Results/", full_scen[[1]], "_dist_par.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[1]] %>%
#   group_by(c(Mother_ID)) %>% # 
#   ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), y = dist_par)) +  
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   #theme_minimal() +  # set theme
#   theme_bw() +  # set theme
#   geom_boxplot(fill="darkolivegreen4") +
#   geom_jitter(color = "grey", fill = "black", width = 0.3) +
#   geom_text(data = . %>% count(Mother_ID), aes(label = n, y = 645), vjust = -0.5) + 
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)")
# 
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[2]], "_dist_par.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[2]] %>%
#   group_by(c(Mother_ID)) %>% # 
#   ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), y = dist_par)) +  
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   #theme_minimal() +  # set theme
#   theme_bw() +  # set theme
#   geom_boxplot(fill="darkolivegreen4") +
#   geom_jitter(color = "grey", fill = "black", width = 0.3) +
#   geom_text(data = . %>% count(Mother_ID), aes(label = n, y = 645), vjust = -0.5) + 
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)")
# 
# dev.off()
# 
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[3]], "_dist_par.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[3]] %>%
#   group_by(c(Mother_ID)) %>% # 
#   ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), y = dist_par)) +  
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   #theme_minimal() +  # set theme
#   theme_bw() +  # set theme
#   geom_boxplot(fill="darkolivegreen4") +
#   geom_jitter(color = "grey", fill = "black", width = 0.3) +
#   geom_text(data = . %>% count(Mother_ID), aes(label = n, y = 645), vjust = -0.5) + 
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)")
# 
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[4]], "_dist_par.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[4]] %>%
#   group_by(c(Mother_ID)) %>% # 
#   ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), y = dist_par)) +  
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   #theme_minimal() +  # set theme
#   theme_bw() +  # set theme
#   geom_boxplot(fill="darkolivegreen4") +
#   geom_jitter(color = "grey", fill = "black", width = 0.3) +
#   geom_text(data = . %>% count(Mother_ID), aes(label = n, y = 645), vjust = -0.5) + 
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)")
# 
# dev.off()

#loop over for all cases 
# 
# for(scen in 1:length(full_scen)){
#   par_scen_df[[scen]] <- par_scen_df[[scen]] %>%
#   mutate(`Parents Are` = case_when(Half_Sibs == FALSE ~ "Not Half Siblings",
#                                                 TRUE ~ "Half Siblings", 
#                                                 NA ~ "No Accession Number")
#   )
# }
# 
# #graph of half-sibling matings group by maternal ID
# png(paste0("Results/Parentage_Figures/", full_scen[[1]], "_halfsib_dist.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[1]] %>%
#   group_by(`Parents Are`) %>%  # group by half siblings to compare the status
#   ggplot(aes(x = Mother_ID, y = dist_par, color = `Parents Are`)) +  
#   geom_jitter(width = 0.2) +
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   scale_color_manual(values = c("cadetblue", "navy")) +
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw()
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[2]], "_halfsib_dist.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[2]] %>%
#   group_by(`Parents Are`) %>%  # group by half siblings to compare the status
#   ggplot(aes(x = Mother_ID, y = dist_par, color = `Parents Are`)) +  
#   geom_jitter(width = 0.2) +
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   scale_color_manual(values = c("cadetblue", "navy")) +
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw()
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[3]], "_halfsib_dist.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[3]] %>%
#   group_by(`Parents Are`) %>%  # group by half siblings to compare the status
#   ggplot(aes(x = Mother_ID, y = dist_par, color = `Parents Are`)) +  
#   geom_jitter(width = 0.2) +
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   scale_color_manual(values = c("cadetblue", "navy")) +
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw()
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[4]], "_halfsib_dist.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[4]] %>%
#   group_by(`Parents Are`) %>%  # group by half siblings to compare the status
#   ggplot(aes(x = Mother_ID, y = dist_par, color = `Parents Are`)) +  
#   geom_jitter(width = 0.2) +
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   scale_color_manual(values = c("cadetblue", "navy")) +
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw()
# dev.off()
# 
# #loop to add hybrid status to data frame 
# for(scen in 1:length(full_scen)){
#   
#   par_scen_df[[scen]] <- par_scen_df[[scen]] %>%
#     mutate(`Offspring Hybrid Status` = case_when(Hybrid_Status == FALSE ~ "Not Hybrid",
#                                                  TRUE ~ "Hybrid"))
#   
# }
# 
# #####barplots of the percentage of hybrids 
# png(paste0("Results/Parentage_Figures/", full_scen[[1]],"_hybrid_boxplot.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[1]] %>%
# ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), 
#                                     y = dist_par, 
#                                   fill = `Offspring Hybrid Status`)) +  
#   geom_boxplot() +
#   scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw() 
# 
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[2]],"_hybrid_boxplot.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[2]] %>%
#   ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), 
#              y = dist_par, 
#              fill = `Offspring Hybrid Status`)) +  
#   geom_boxplot() +
#   scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw() 
# 
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[3]],"_hybrid_boxplot.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[3]] %>%
#   ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), 
#              y = dist_par, 
#              fill = `Offspring Hybrid Status`)) +  
#   geom_boxplot() +
#   scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw() 
# 
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[4]],"_hybrid_boxplot.png"), 
#     res = 600, width = 5000, height = 3500)
# par_scen_df[[4]] %>%
#   ggplot(aes(x = fct_rev(fct_infreq(Mother_ID)), 
#              y = dist_par, 
#              fill = `Offspring Hybrid Status`)) +  
#   geom_boxplot() +
#   scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
#   expand_limits(y = c(0, 650)) +  # set limits for graph
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw() 
# 
# dev.off()
# 
# 
# ##create a list for 
# #table to present the candidate fathers 
# species_count_list <- list()
# 
# #loop to create species count lists for parent assignments
# for(scen in 1:length(full_scen)){
#   
#   #create a data frame 
#   species_count_list[[scen]] <- as.data.frame(table(par_scen_df[[scen]]$Candidate_Father_Species))
#   
#   #rename columns 
#   names(species_count_list[[scen]]) <- c("Species", "Count")
#   
#   #organize data frame 
#   species_count_list[[scen]]$Species <- factor(species_count_list[[scen]]$Species, 
#                                                levels=species_count_list[[scen]]$Species[order(-species_count_list[[scen]]$Count)])
#   
#   
# }
# 
# 
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[2]], "_hybrid_barplot.png"),
#     res = 600, width = 5000, height = 3500)
# species_count_list[[2]] %>%
#   ggplot(aes(x = Species, y=Count))+
#   geom_bar(stat = "identity", fill = "darkgreen") + 
#   labs(title="Count of Candidate Father Trees per Species", 
#        x="Candidate Father Species") +
#   theme_bw()
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[3]], "_hybrid_barplot.png"),
#     res = 600, width = 5000, height = 3500)
# species_count_list[[3]] %>%
#   ggplot(aes(x = Species, y=Count))+
#   geom_bar(stat = "identity", fill = "darkgreen") + 
#   labs(title="Count of Candidate Father Trees per Species", 
#        x="Candidate Father Species") +
#   theme_bw()
# dev.off()
# 
# png(paste0("Results/Parentage_Figures/", full_scen[[4]], "_hybrid_barplot.png"),
#     res = 600, width = 5000, height = 3500)
# species_count_list[[4]] %>%
#   ggplot(aes(x = Species, y=Count))+
#   geom_bar(stat = "identity", fill = "darkgreen") + 
#   labs(title="Count of Candidate Father Trees per Species", 
#        x="Candidate Father Species") +
#   theme_bw()
# dev.off()
# 
# ################## Unused Loops ------------- 
# for(scen in seq_along(full_scen)){
#   
#   mat_off <- par_scen_df[[scen]] %>%
#     ggplot() +
#     geom_bar(aes(y = sort(Candidate_father_ID))) +
#     facet_wrap(~`Mother_ID`) + 
#     scale_x_continuous(n.breaks = 9) +
#     labs(title = "Count of Offspring per Candidate Father and Maternal Tree Pairs", 
#          y = "Candidate Father ID", x = "Count of Offspring")
#   
#   png(paste0("Results/Parentage_Figures/", full_scen[[scen]], "_mat_off.png"), 
#       width = 5000, height = 3500, res = 600)
#   
#   mat_off
#   
# }
# dev.off()
# 

