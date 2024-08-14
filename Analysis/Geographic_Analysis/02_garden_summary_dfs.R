#####################
#     Libraries     #
#####################


###########################
#     Load Data Files     #
###########################
#list out garden data files 
setwd("../..")

#list out garden occ files 
garden_list_occ <- list.files(pattern = "_occ_df.csv$",
                              path = "Results/Geographic_Analyses")

#######################################
#     Create Summary Data Frames      #
#######################################

for(garden in seq_along(garden_list_occ)){
  
  sp_df <- read.csv(paste0("./Results/Geographic_Analyses/",
                           garden_list_occ[[garden]]))
  
  species_list <- unique(sp_df$Species_Name)
  
  garden_species_df <- matrix(nrow = length(species_list), ncol = 1)
  
  rownames(garden_species_df) <- species_list
  
  for(sp in 1:length(species_list)){
    
    garden_species_df[sp,1] <- length(sp_df[sp_df$Species_Name == species_list[[sp]],][,1])
   
    write.csv(garden_species_df, 
              paste0("./Results/Geographic_Analyses/",gsub("\\..*",'', garden_list_occ[[garden]]), "_species_df.csv"))
    
  }
  
}
