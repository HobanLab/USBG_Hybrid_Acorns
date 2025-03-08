###This script creates summarization files and figures for the final publication
##of this study. Summary data file code was created by Emily Schumacher and 
#figure code was created by Mikaely Evans. 

#####################
#     Libraries     #
#####################

library(tidyverse)

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
par_sum_df <- list(par_sum_df[[3]], par_sum_df[[1]],
                   par_sum_df[[4]], par_sum_df[[2]])

#list out analysis data frames 
par_sum_analysis <- list.files(path = "Results/Parentage_Results/CSV_Files",
                               pattern = "analysis_df")

#order list 
par_sum_analysis <- c(par_sum_analysis[[1]], par_sum_analysis[[2]],
                      par_sum_analysis[[4]], par_sum_analysis[[3]])

#load full scenario data frames 
full_scen <- c("all_loci_AF", #all loci included with all father assignments 
               "all_loci_HCF", #all loci with only high confidence fathers included
               "red_loci_AF", #reduced loci with all father assignments
               "red_loci_HCF" #reduce loci with only high confidence father assignments included
)

#pull in UHA database 
UHA_database <- read.csv("./Data_Files/CSV_Files/UHA_database_clean.csv")

#save out list of maternal IDs 
mat_n_ids <- paste0("MT",1:11)

#list out the mt ids 
mat_ids <- list()

for(mat in mat_n_ids){
  
  mat_ids[[mat]] <- UHA_database[UHA_database$MT_ID == mat,"Tissue_ID"]
  
}

#create MT 
UHA_res_df$MT_ID <- dplyr::case_when(
  UHA_res_df$Maternal_ID == "UHA-0009" ~ "MT1",
  UHA_res_df$Maternal_ID == "UHA-0010" ~ "MT2",
  UHA_res_df$Maternal_ID == "UHA-0012" ~ "MT3",
  UHA_res_df$Maternal_ID == "UHA-0013" ~ "MT4",
  UHA_res_df$Maternal_ID == "UHA-0014" ~ "MT5",
  UHA_res_df$Maternal_ID == "UHA-0015" ~ "MT6",
  UHA_res_df$Maternal_ID == "UHA-0016" ~ "MT7",
  UHA_res_df$Maternal_ID == "UHA-0257" ~ "MT8",
  UHA_res_df$Maternal_ID == "UHA-0260" ~ "MT9",
  UHA_res_df$Maternal_ID == "UHA-0261" ~ "MT10",
  UHA_res_df$Maternal_ID == "UHA-0351" ~ "MT11"
)

#save as a factor 
UHA_res_df$MT_ID <- as.factor(UHA_res_df$MT_ID)
UHA_res_df$MT_ID <- factor(UHA_res_df$MT_ID, 
                           levels=c("MT1", "MT2", "MT3", "MT4",
                                    "MT5","MT6","MT7","MT8",
                                    "MT9","MT10","MT11"))

## add candidate father names 
UHA_res_df$CF_ID <- dplyr::case_when(
  UHA_res_df$Candidate_father_ID == "UHA-0011" ~ "CF1",
  UHA_res_df$Candidate_father_ID == "UHA-0013" ~ "CF2",
  UHA_res_df$Candidate_father_ID == "UHA-0015" ~ "CF3",
  UHA_res_df$Candidate_father_ID == "UHA-0177" ~ "CF4",
  UHA_res_df$Candidate_father_ID == "UHA-0224" ~ "CF5",
  UHA_res_df$Candidate_father_ID == "UHA-0256" ~ "CF6",
  UHA_res_df$Candidate_father_ID == "UHA-0257" ~ "MT8",
  UHA_res_df$Candidate_father_ID == "UHA-0258" ~ "CF7",
  UHA_res_df$Candidate_father_ID == "UHA-0259" ~ "CF8",
  UHA_res_df$Candidate_father_ID == "UHA-0260" ~ "MT9",
  UHA_res_df$Candidate_father_ID == "UHA-0261" ~ "MT10",
  UHA_res_df$Candidate_father_ID == "UHA-0262" ~ "CF9",
  UHA_res_df$Candidate_father_ID == "UHA-0263" ~ "CF10",
  UHA_res_df$Candidate_father_ID == "UHA-0264" ~ "CF11",
  UHA_res_df$Candidate_father_ID == "UHA-0265" ~ "CF12",
  UHA_res_df$Candidate_father_ID == "UHA-0266" ~ "CF13",
  UHA_res_df$Candidate_father_ID == "UHA-0267" ~ "CF14",
  UHA_res_df$Candidate_father_ID == "UHA-0268" ~ "CF15",
  UHA_res_df$Candidate_father_ID == "UHA-0269" ~ "CF16",
  UHA_res_df$Candidate_father_ID == "UHA-0293" ~ "CF17",
  UHA_res_df$Candidate_father_ID == "UHA-0303" ~ "CF18",
  UHA_res_df$Candidate_father_ID == "UHA-0351" ~ "MT11",
  UHA_res_df$Candidate_father_ID == "UHA-0476" ~ "CF19",
  UHA_res_df$Candidate_father_ID == "UHA-0512" ~ "CF20"
  
  
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
                                     par_sum_df[[par]]))[,-1]
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

## create a summary data frame with the maternal IDs
par_sum_stat_df <- matrix(nrow = length(mat_ids),
                          ncol = 11)

#add the maternal ids in a column
par_sum_stat_df[,1] <-  as.character(mat_ids)

#loop over all parentage scenarios with a summary data file 
for(par in seq_along(par_sum_df)){
  
  temp_df <- read.csv(paste0("./Results/Parentage_Results/CSV_Files/",par_sum_analysis[[par]]))
  
  #loop over each maternal individual to create data frame 
  for(mat in seq_along(mat_ids)){
    
    #add a column for the offspring counts by mom 
    par_sum_stat_df[mat,2] <- length(temp_df[temp_df$Maternal_ID == mat_ids[[mat]],"Offspring_ID"])
    
    #add a column for number of dads per mom 
    par_sum_stat_df[mat,3] <- length(unique(temp_df[temp_df$Maternal_ID == mat_ids[[mat]],"Candidate_father_ID"]))
    
    #add columns with the number of hybrid per mom and the number of hybrid fathers
    par_sum_stat_df[mat,4] <- length(temp_df[(temp_df$Maternal_ID == mat_ids[[mat]]) & (temp_df$Hybrid_Status == TRUE),"Candidate_father_ID"])
    par_sum_stat_df[mat,5] <- length(unique(temp_df[(temp_df$Maternal_ID == mat_ids[[mat]]) & (temp_df$Hybrid_Status == TRUE),"Candidate_father_ID"]))
  
    #add half-sib columns 
    par_sum_stat_df[mat,6] <- length(temp_df[(temp_df$Maternal_ID == mat_ids[[mat]]) & (temp_df$Half_Sibs == TRUE),"Candidate_father_ID"])
    par_sum_stat_df[mat,7] <- length(unique(temp_df[(temp_df$MT_ID == mat_ids[[mat]]) & (temp_df$Half_Sibs == TRUE),"Candidate_father_ID"]))
    
      if(par_sum_stat_df[mat,2] == 0){
        par_sum_stat_df[mat,8] <- 0
        par_sum_stat_df[mat,9] <- 0
        par_sum_stat_df[mat,10] <- 0
        par_sum_stat_df[mat,11] <- 0
        
      }else{
      
        par_sum_stat_df[mat,8] <- min(temp_df[(temp_df$Maternal_ID == mat_ids[[mat]]),"dist_par"])
        par_sum_stat_df[mat,9] <- max(temp_df[(temp_df$Maternal_ID == mat_ids[[mat]]),"dist_par"])
        
        #recode hybrid status of min distance 
        par_sum_stat_df[mat,10] <- length(temp_df[temp_df$dist_par == min(temp_df[(temp_df$Maternal_ID == mat_ids[[mat]]),"dist_par"]) & temp_df$Hybrid_Status == TRUE,"Candidate_Father_Species"])
        #recode hybrid status of min distance 
        par_sum_stat_df[mat,11] <- length(temp_df[temp_df$dist_par == max(temp_df[(temp_df$Maternal_ID == mat_ids[[mat]]),"dist_par"]) & temp_df$Hybrid_Status == TRUE,"Candidate_Father_Species"])
        
    }
     
    # #add colnames
    colnames(par_sum_stat_df) <- c("MT_ID", "Off_N",
                                  "Fathers_N", "Hybrid_Off_N",
                                  "Hybrid_Father_N", "Half_Sib_Off_N",
                                  "Half_Sib_Father_N", "Min_Dist", "Max_Dist",
                                  "Hybrid_Min", "Hybrid_Max")
    # #write out summary data frame with scenario 
    write.csv(par_sum_stat_df, paste0("Results/Parentage_Results/CSV_Files/", full_scen[[par]], '_par_sum_stat_df.csv'),
             row.names = FALSE)
  }
  
  
  
}

###################
#     Figures     #
###################

###### Figure 1 - Species Count Barplot -------------------

#create summary father data frame 
UHA_father_df <- as.data.frame(table(UHA_res_df$Candidate_Father_Species))
#rename 
colnames(UHA_father_df) <- c("Species", "Count")

#resort data frame 
UHA_father_df <- UHA_father_df %>%
                   dplyr::arrange(across("Count",desc))

#### Now write out figures 
pdf(paste0("Results/Parentage_Results/Figures/AL_HCF_hybrid_species_count.pdf"),
    height = 6, width = 8)

UHA_father_df %>%
  ggplot(aes(x = fct_rev(fct_reorder(Species, Count)), y = Count)) +
  geom_col(fill="forestgreen") +
  geom_text(aes(label = paste("n =", Count)), vjust = -0.5) + 
  scale_y_continuous(limits = c(0,150)) +
  theme_bw() +
  labs(x="Candidate Father Species",y = "Number of Offspring") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 1))

dev.off()

###### Figure 2 - Maternal/Paternal Tree Distances Boxplot -------------------
UHA_res_df$hybrid_update <- NA
UHA_res_df[UHA_res_df$Hybrid_Status == "TRUE","Hybrid Status"] <- "Hybrid"
UHA_res_df[UHA_res_df$Hybrid_Status == "FALSE","Hybrid Status"] <- "Not Hybrid"

## write out figure text to compare hybrid status/distance by mother
pdf(paste0("Results/Parentage_Results/Figures/AL_HCF_dist_par_hybrid.pdf"),
    height = 6, width = 8)

UHA_res_df %>%
  ggplot(aes(x = MT_ID, 
             y = dist_par)) +  
  geom_boxplot(fill = "darkolivegreen4") +
  geom_jitter(aes(fill = `Hybrid Status`), width = 0.2, size = 3, shape = 21, color = "black") +
  scale_y_continuous(limits = c(0,1000)) +  # set limits for graph
  scale_fill_manual(values = c("deeppink", "grey")) +
  geom_text(data = . %>% count(MT_ID), 
            aes(label = paste("n =", n), y = 665), vjust = -3) + 
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5))


dev.off()

###### Figure 3 - Boxplot Relatedness Parent Boxplot ------------------

#replace column text for halfsibs
UHA_res_df[["Parents Are"]] <- dplyr::case_when(
  UHA_res_df$Half_Sibs == TRUE ~ "Half Siblings",
  UHA_res_df$Half_Sibs == FALSE ~ "Not Half Siblings"
)
UHA_res_df[is.na(UHA_res_df$`Parents Are`),]$`Parents Are` <- "Unknown"

#graph of half-sibling matings group by maternal ID
pdf(paste0("Results/Parentage_Results/Figures/AL_HCF_dist_par_halfsib_relanl.pdf"),
    width = 12, height = 8)
UHA_res_df %>%
  group_by(sum_relate) %>%  # group by half siblings to compare the status
  ggplot(aes(x = MT_ID, y = dist_par, fill = sum_relate)) +  
  geom_boxplot(fill = "gray92")+
  geom_jitter(aes(fill = sum_relate), width = 0.2, size = 5, shape = 21, color = "black") +
  expand_limits(y = c(0, 800)) +  # set limits for graph
  scale_fill_manual(values = c("moccasin", "cornflowerblue","white")) +
  xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.3))
dev.off()


#### Supplemental Figures - Not in manuscript --------------

##Barplot breakdown of offspring count 

#replace column text for candidate father species
UHA_res_df[["Candidate Father Species"]] <- as.factor(UHA_res_df$Candidate_Father_Species)

#make a factor column to enable easier counting of the fathers
UHA_res_df$CF_ID <- as.factor(UHA_res_df$CF_ID)

#write out count by father figure
png(paste0("Results/Parentage_Results/Figures/AL_HCF_count_by_father.png"),
    res = 600, width = 7500, height = 3500)

UHA_res_df %>%
  ggplot(aes(x = fct_infreq(CF_ID))) +
  geom_bar(aes(fill = `Candidate Father Species`)) +
  scale_fill_manual(values=c("deepskyblue3", "goldenrod3","darkseagreen4","lightgoldenrod1")) +
  geom_text(data = . %>% count(CF_ID), 
            aes(label = paste("n =", n), y = 31), vjust = -3,
            size = 3) +
  scale_y_continuous(limits = c(0,35)) +  # set limits for graph
  xlab("Candidate Father ID") + ylab("Offspring Count") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.3))
dev.off()

### Barplots of assigned fathers by mother 
png("Results/Parentage_Results/Figures/AL_HCF_mat_count_CF.png",
    res = 600, width = 5000, height = 6500)
UHA_res_df %>%
  ggplot() +
  geom_bar(aes(y = sort(CF_ID))) +
  facet_wrap(~`MT_ID`) +
  scale_x_continuous(n.breaks = 9) +
  labs(y = "Candidate Father ID", x = "Count of Offspring") +
  theme_bw()
dev.off()

########### Make DBH figures

#create count df
cf_list <- unique(UHA_res_df$CF_ID)

#set up data frame
count_df <- as.data.frame(matrix(nrow = length(cf_list),
                                 ncol = 4))
#first column - candidate father list
count_df[,1] <- cf_list

#loop over all cases 
for(cf in seq_along(cf_list)){
  
  #second column - count 
  count_df[cf,2] <- length(UHA_res_df[UHA_res_df$CF_ID == cf_list[[cf]],"Offspring_ID"])
  
  #add in dbh
  count_df[cf,3] <- unique(UHA_res_df[UHA_res_df$CF_ID == cf_list[[cf]],"DBH_avg"])
  
  #add in species ID
  count_df[cf,4] <- unique(UHA_res_df[UHA_res_df$CF_ID == cf_list[[cf]],"Candidate_Father_Species"])
  
}

#name columns 
colnames(count_df) <- c("Candidate_Father_ID", "Count", "DBH_avg", "Species")
#coerce to data frame
count_df <- as.data.frame(count_df)
class(count_df$DBH_avg) <- "numeric"

#now sort by dbh
count_df <- count_df %>%
              dplyr::arrange(across("DBH_avg",desc))
count_df <- as.data.frame(count_df)


### DBH Figure 

#plot linear regression - loop over numeric
species_list <- c("Quercus muehlenbergii","Quercus alba","Quercus macrocarpa","Quercus prinoides")
#create color list 
color_list <- c("darkseagreen4", "deepskyblue3", "goldenrod3", "lightgoldenrod1")

#count df color 
count_df$color <- NA

for(sp in seq_along(species_list)){
  
  count_df[count_df$Species == species_list[[sp]],"color"] <- color_list[[sp]]
  
}

#save out linear regression info 
reg <- lm(count_df$Count ~ count_df$DBH_avg)
sum_reg <- summary(reg)
#save out r2 and p value 
adj_r2 <- sum_reg$adj.r.squared
pvalue <- sum_reg$coefficients[[8]]

#create legend
rp <- vector('expression',2)
rp[1] <- substitute(expression(italic(R)^2 == MYVALUE), 
                            list(MYVALUE = format(adj_r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                           list(MYOTHERVALUE = format(pvalue, digits = 2)))[2]
#write out 
png("Results/Parentage_Results/Figures/AL_HCF_dbh_offspring_count.png",
    res = 600, width = 5000, height = 3500)
#plot linear regression based on color
plot(x = count_df[count_df$Species == species_list[[1]],]$DBH_avg, 
     y = count_df[count_df$Species == species_list[[1]],]$Count,
     pch = 17, 
     col = count_df[count_df$Species == species_list[[1]],]$color,
     ylim = c(0,30),
     xlab = "Average DBH (cm)",
     ylab = "Offspring Count",
     xlim = c(0,100))
points(x = count_df[count_df$Species == species_list[[2]],]$DBH_avg, 
       y = count_df[count_df$Species == species_list[[2]],]$Count,
       pch = 17, 
       col = count_df[count_df$Species == species_list[[2]],]$color)
points(x = count_df[count_df$Species == species_list[[3]],]$DBH_avg, 
       y = count_df[count_df$Species == species_list[[3]],]$Count,
       pch = 17, 
       col = count_df[count_df$Species == species_list[[3]],]$color)
points(x = count_df[count_df$Species == species_list[[4]],]$DBH_avg, 
       y = count_df[count_df$Species == species_list[[4]],]$Count,
       pch = 17, 
       col = count_df[count_df$Species == species_list[[4]],]$color)

#add regression 
abline(reg, col = "blue", lwd = 2)

#add legend - r2 and pvalue
legend("topright", legend = rp, bty = 'n', border = "black", 
       pt.cex = 1.5, cex = 0.8, pch = 15, col = "blue")
#add species legend
legend('topleft', legend = species_list,
       pt.cex = 1.5, cex = 0.8,
       pch = 17, col = color_list)
dev.off()

############## Not currently in use ----------------------

##Histogram of DBHs by species 
#create a plot of the DBH across fathers 
png("Results/Parentage_Results/Figures/AL_HCF_barplot_DBH.png",
    res = 600, width = 5200, height = 3500)

count_df %>%
  ggplot(aes(x = Candidate_Father_ID, y = DBH_avg, fill = Species)) +
  scale_fill_manual(values=c("deepskyblue3", "goldenrod3","darkseagreen4","lightgoldenrod1")) +
  geom_bar(stat = "identity") +
  xlab("Candidate Father ID") + ylab("DBH (cm)") +
  theme_bw() +
  scale_y_continuous(limits = c(0,150))

dev.off()


# 
# #graph of half-sibling matings group by maternal ID
# png(paste0("Results/Parentage_Results/Figures/AL_HCF_dist_par_halfsib.png"),
#     res = 600, width = 5200, height = 3500)
# UHA_res_df %>%
#   group_by(`Parents Are`) %>%  # group by half siblings to compare the status
#   ggplot(aes(x = MT_ID, y = dist_par, fill = `Parents Are`)) +  
#   geom_boxplot(fill = "azure2")+
#   geom_jitter(aes(fill = `Parents Are`), width = 0.2, size = 3, shape = 21, color = "black") +
#   expand_limits(y = c(0, 800)) +  # set limits for graph
#   scale_fill_manual(values = c("cadetblue", "navy","white")) +
#   xlab("Maternal Tree ID") + ylab("Distance between parents (m)") +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
#         axis.text.y = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16),
#         plot.title = element_text(size = 18, hjust = 0.3))
# dev.off()
# ###### Boxplot of distances between parents 
# png(paste0("Results/Parentage_Results/", full_scen[[1]], "_dist_par.png"), 
#     res = 600, width = 5000, height = 3500)
# par_sum_df[[1]] %>%
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
# par_sum_df[[2]] %>%
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
# par_sum_df[[3]] %>%
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
# par_sum_df[[4]] %>%
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
#   par_sum_df[[scen]] <- par_sum_df[[scen]] %>%
#   mutate(`Parents Are` = case_when(Half_Sibs == FALSE ~ "Not Half Siblings",
#                                                 TRUE ~ "Half Siblings", 
#                                                 NA ~ "No Accession Number")
#   )
# }
# 
# #graph of half-sibling matings group by maternal ID
# png(paste0("Results/Parentage_Figures/", full_scen[[1]], "_halfsib_dist.png"), 
#     res = 600, width = 5000, height = 3500)
# par_sum_df[[1]] %>%
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
# par_sum_df[[2]] %>%
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
# par_sum_df[[3]] %>%
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
# par_sum_df[[4]] %>%
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
#   par_sum_df[[scen]] <- par_sum_df[[scen]] %>%
#     mutate(`Offspring Hybrid Status` = case_when(Hybrid_Status == FALSE ~ "Not Hybrid",
#                                                  TRUE ~ "Hybrid"))
#   
# }
# 
# #####barplots of the percentage of hybrids 
# png(paste0("Results/Parentage_Figures/", full_scen[[1]],"_hybrid_boxplot.png"), 
#     res = 600, width = 5000, height = 3500)
# par_sum_df[[1]] %>%
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
# par_sum_df[[2]] %>%
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
# par_sum_df[[3]] %>%
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
# par_sum_df[[4]] %>%
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
#   species_count_list[[scen]] <- as.data.frame(table(par_sum_df[[scen]]$Candidate_Father_Species))
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
#   mat_off <- 
#   png(paste0("Results/Parentage_Figures/", full_scen[[scen]], "_mat_off.png"), 
#       width = 5000, height = 3500, res = 600)
#   
#   mat_off
#   
# }
# dev.off()
# 

