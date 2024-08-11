#####################
#     Libraries     #
#####################

library(tidyverse)
library(sjmisc)

###########################
#     Load Data Files     #
###########################
#setwd("USBG_Hybrid_Acorn/Data_Files/Geographic_Files")

#
QUoc <- list.files(pattern = "_QUocc.csv")

red_oak_names <- c("Quercus acerifolia", "Quercus acutifolia", "Quercus agrifolia", 
                   "Quercus albocincta", "Quercus aristata", "Quercus arkansana", 
                   "Quercus buckleyi", "Quercus canbyi", "Quercus candicans", 
                   "Quercus castanea", "Quercus coccinea", "Quercus coccolobifolia", 
                   "Quercus coffeicolor", "Quercus conspersa", "Quercus costaricensis", 
                   "Quercus crassifolia", "Quercus crassipes", "Quercus cualensis", 
                   "Quercus delgadoana", "Quercus depressa", "Quercus durifolia", 
                   "Quercus dysophylla", "Quercus eduardii", "Quercus ellipsoidalis", 
                   "Quercus elliptica", "Quercus emoryi", "Quercus falcata", 
                   "Quercus frutex", "Quercus gentryi", "Quercus glabrescens", 
                   "Quercus gravesii", "Quercus graciliformis", "Quercus georgiana", 
                   "Quercus hemisphaerica", "Quercus hintonii", "Quercus hintoniorum", 
                   "Quercus hirtifolia", "Quercus humboldtii", "Quercus hypoleucoides", 
                   "Quercus hypoxantha", "Quercus ilicifolia", "Quercus iltisii", 
                   "Quercus imbricaria", "Quercus incana", "Quercus inopina", 
                   "Quercus jonesii", "Quercus kelloggii", "Quercus laevis", 
                   "Quercus laurifolia", "Quercus laurina", "Quercus marilandica", 
                   "Quercus martinezii", "Quercus mcvaughii", "Quercus mexicana", 
                   "Quercus myrtifolia", "Quercus nigra", "Quercus pagoda", 
                   "Quercus parvula", "Quercus palustris", "Quercus phellos",
                   "Quercus planipocula", "Quercus potosina", "Quercus pumila", 
                   "Quercus radiata", "Quercus rapurahuensis", "Quercus resinosa", 
                   "Quercus robusta", "Quercus rysophylla", "Quercus rubra", 
                   "Quercus salicifolia", "Quercus sapotifolia", "Quercus scytophylla", 
                   "Quercus shumardii", "Quercus sideroxyla", "Quercus splendens", 
                   "Quercus skinneri", "Quercus tardifolia", "Quercus texana", 
                   "Quercus urbanii", "Quercus uxoris", "Quercus velutina", 
                   "Quercus viminea", "Quercus wislizeni", "Quercus xalapensis")

white_oak_names <- c("Quercus ajoensis", "Quercus alba", "Quercus aliena", 
                     "Quercus arizonica", "Quercus austrina", "Quercus ? bebbiana", 
                     "Quercus berberidifolia", "Quercus bicolor", "Quercus ? bimundorum", 
                     "Quercus boyntonii", "Quercus brandegeei", "Quercus carmenensis", 
                     "Quercus chapmanii", "Quercus chihuahuensis", 
                     "Quercus cornelius-mulleri", "Quercus conzattii", "Quercus copeyensis", 
                     "Quercus dalechampii", "Quercus depressipes", "Quercus deserticola", 
                     "Quercus diversifolia", "Quercus douglasii", "Quercus dumosa",
                     "Quercus durata", "Quercus engelmannii", "Quercus fabrei", 
                     "Quercus faginea", "Quercus fulva", "Quercus furuhjelmi", 
                     "Quercus fusiformis", "Quercus gambelii", "Quercus garryana", 
                     "Quercus geminata", "Quercus glaucoides", "Quercus greggii", 
                     "Quercus griffithii", "Quercus grisea", "Quercus hartwissiana", 
                     "Quercus havardii", "Quercus hiholensis", "Quercus hinckleyi", 
                     "Quercus hondurensis", "Quercus infectoria", "Quercus infectoria", 
                     "Quercus insignis", "Quercus intricata", "Quercus john-tuckeri", 
                     "Quercus laceyi", "Quercus laeta", "Quercus lanata", 
                     "Quercus lancifolia", "Quercus leucotrichophora", "Quercus liebmannii",
                     "Quercus lobata", "Quercus lusitanica", "Quercus lyrata", 
                     "Quercus macrocarpa", "Quercus margarettae", "Quercus magnoliifolia", 
                     "Quercus martinezii", "Quercus mohriana", "Quercus montana", 
                     "Quercus michauxii", "Quercus microphylla", "Quercus minima", 
                     "Quercus mongolica", "Quercus muehlenbergii", "Quercus oblongifolia", 
                     "Quercus obtusata", "Quercus oglethorpensis", "Quercus oleoides", 
                     "Quercus oocarpa", "Quercus pacifica", "Quercus peduncularis", 
                     "Quercus petraea", "Quercus polymorpha", "Quercus praeco", 
                     "Quercus prinoides", "Quercus pubescens", "Quercus pungens", 
                     "Quercus robur", "Quercus rugosa", "Quercus sagraeana", 
                     "Quercus ? schuettei", "Quercus sebifera", "Quercus serrata", 
                     "Quercus similis", "Quercus sinuata", "Quercus sinuata var. breviloba", 
                     "Quercus sinuata var. sinuata", "Quercus stellata", "Quercus striatula",
                     "Quercus subspathulata", "Quercus tarahumara", "Quercus toumeyi", 
                     "Quercus tuberculata", "Quercus turbinella", "Quercus ? turneri", 
                     "Quercus vaseyana", "Quercus vincentensis", "Quercus virginiana", 
                     "Quercus welshii", "Quercus wutaishanica")


intermediate_oaks <- c("Quercus cedrosensis", "Quercus chrysolepis", "Quercus palmeri",
                       "Quercus tomentella", "Quercus vacciniifolia")

cerris_oaks <- c("Quercus acutissima", "Quercus alnifolia", "Quercus aquifolioides",
                 "Quercus brandisiana", "Quercus brantii", "Quercus calliprinos",
                 "Quercus castaneifolia", "Quercus cerris", "Quercus chenii", 
                 "Quercus coccifera", "Quercus floribunda", "Quercus franchetii",
                 "Quercus ilex", "Quercus ithaburensis", "Quercus libani",
                 "Quercus macrolepis", "Quercus miyagii", "Quercus pannosa",
                 "Quercus semecarpifolia", "Quercus spinosa", "Quercus suber",
                 "Quercus trojana", "Quercus variabilis")

###################################
#     Organizing Data Frames      #
###################################
#create list to store duplicates 
sp_limited_names <- list()

garden <- 3

for(garden in 1:length(QUoc)){
  
  #load in each garden data frame 
  garden_df <- read.csv(QUoc[[garden]])
  
  #create a column for red oaks 
  garden_df$red_oak_t <- str_contains(red_oak_names, garden_df$Species_Name)
  
  #create a column for white oaks 
  garden_df$white_oak_t <- str_contains(white_oak_names, garden_df$Species_Name)
  
  #create a column for intermediate oaks 
  garden_df$intermediate_oaks_t <- str_contains(intermediate_oaks, garden_df$Species_Name)
  
  #create a column for cerris oaks 
  garden_df$cerris_oaks_t <- str_contains(cerris_oaks, garden_df$Species_Name)
  
  ##now create a section column 
  garden_df$Section <- NA
  
  #input section names 
  garden_df[garden_df$red_oak_t == TRUE,]$Section <- "red_oak"
  garden_df[garden_df$white_oak_t == TRUE,]$Section <- "white_oak"
  garden_df[garden_df$intermediate_oaks_t == TRUE,]$Section <- "intermediate_oak"
  garden_df[garden_df$cerris_oaks_t == TRUE,]$Section <- "cerris_oak"
  
  #write out data frames with the sections added
  write.csv(garden_df, paste0("../../Analyses/", gsub("_.*",'', QUoc[[garden]]),
                              "_occ_df.csv"))
  
  dup_df <- garden_df[duplicated(garden_df$Species_Name) == TRUE,]
  
  #create a list to save duplicates 
  dup_sp_list <- unique(dup_df$Species_Name)
  
  garden_dup_limit_df <- matrix(nrow = length(dup_sp_list), ncol = 1)
  
  #loop to check for duplicates  
  for(dup in 1:length(dup_sp_list)){
     
    #now limit the garden df by the species names 
    garden_dup_df <- garden_df[garden_df$Species_Name == dup_sp_list[[dup]],]
    
    garden_dup_limit_df[dup,1] <- length(garden_dup_df[,1])
    
    rownames(garden_dup_limit_df) <- dup_sp_list
  
    sp_limited_names[[garden]] <- names(which(garden_dup_limit_df[,1] >= 10))
    
    #limit occ records by species names 
    garden_df$limit_names <- str_contains(sp_limited_names[[garden]], garden_df$Species_Name)
    
    #limit data frame by this column 
    garden_occ_limited_df <- garden_df[garden_df$limit_names == TRUE,]
    
    write.csv(garden_occ_limited_df, paste0("../../Analyses/Geographic_Results/",
                                            gsub("_.*",'', QUoc[[garden]]), 
                                           "_limited_occ_df.csv"), 
              row.names = FALSE)
      
  }
}







