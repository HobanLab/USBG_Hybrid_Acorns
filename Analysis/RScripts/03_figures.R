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

#load full parentage file
full_parentage <- read_csv("Data_Files/CSV_Files/UHA_full_parentage.csv")

###########################
#     Visualizations      #
###########################

#plot bar graphs of the offspring count for each mother and which individual is the father
full_parentage %>%
  drop_na(`Candidate_father_ID`) %>%
  ggplot() +
  geom_bar(aes(y = `Candidate_father_ID`)) +
  facet_wrap(~`Mother_ID`) + 
  scale_x_continuous(n.breaks = 9) +
  labs(title = "Count of Offspring per Candidate Father and Maternal Tree Pairs", 
       y = "Candidate Father ID", x = "Count of Offspring")


mf_boxplot <- full_parentage %>%
  group_by(c(`Mother_ID`)) %>% # Grouped by mother ID for analysis across maternal tree
  filter(!is.na(`Candidate_father_ID`)) %>%
  ggplot(aes(x = fct_rev(fct_infreq(`Mother_ID`)), y = distance_between_parents)) +  #This is making me think that I actually grouped it by count of occurrences lowest to highest instead of highest average distance between parents. As we discussed with Sean, the order is not that important.
  expand_limits(y = c(0, 650)) +  # set limits for graph
  #theme_minimal() +  # set theme
  theme_bw() +  # set theme
  geom_boxplot(fill="darkolivegreen4") +
  geom_jitter(color = "grey", fill = "black", width = 0.3) +
  geom_text(data = . %>% count(`Mother_ID`), aes(label = n, y = 645), vjust = -0.5) + 
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)")

#boxplot of the mothers and fathers 
png("Results/mothers_fathers.png", res = 600,
    width = 3500, height = 3500)

mf_boxplot

dev.off()

#create a column for the renamed cases for half siblings
full_parentage <- full_parentage %>%
  mutate(`Parents Are` = case_when(Half_Sibs == FALSE ~ "Not Half Siblings",
                     TRUE ~ "Half Siblings")
  )

#graph of half-sibling matings group by maternal ID
png("Results/half_sibs_jitter.png", res = 600,
    width = 5000, height = 3500)
full_parentage %>%
  group_by(`Parents Are`) %>%  # group by half siblings to compare the status
  filter(!is.na(`Candidate_father_ID`)) %>%
  filter(!is.na(`Parents Are`)) %>%
  ggplot(aes(x = Mother_ID, y = distance_between_parents, color = `Parents Are`)) +  
  geom_jitter(width = 0.2) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  scale_color_manual(values = c("cadetblue", "navy")) +
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw()
dev.off()

#organize data frame for hybrid status
full_parentage <- full_parentage %>%
  mutate(`Offspring Hybrid Status` = case_when(Hybrid_Status == FALSE ~ "Not Hybrid",
                                   TRUE ~ "Hybrid"))

#barplot of the percentage of hybrids 
png("Results/hybrid_offpring.png", res = 600,
    width = 5000, height = 3500)
ggplot(data = full_parentage, aes(x = fct_rev(fct_infreq(Mother_ID)), 
                                    y = distance_between_parents, 
                                  fill = `Offspring Hybrid Status`)) +  
  geom_boxplot() +
  scale_fill_manual(values = c("darkseagreen", "darkgreen")) +
  expand_limits(y = c(0, 650)) +  # set limits for graph
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw() 
dev.off()

#table to present the candidate fathers 
species_count <- as.data.frame(table(full_parentage$Candidate_Father_Species))
names(species_count) <- c("Species", "Count")
species_count$Species <- factor(species_count$Species, 
                                levels=species_count$Species[order(-species_count$Count)])
#now plot species count table
species_count %>%
  ggplot(aes(x = Species, y=Count))+
  geom_bar(stat = "identity", fill = "darkgreen") + 
  labs(title="Count of Candidate Father Trees per Species", 
       x="Candidate Father Species") +
  theme_bw()
