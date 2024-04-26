## Project Description
This github project analyzes oak seedlings to detect hybridization using DNA markers. Oak trees are critical ecosystem creators and are of very high conservation value. Additionally, oaks are an exceptional species, meaning they cannot be seedbanked using traditional methods and can only be conserved in living collections like botanic gardens or arboreta. Maintaining living collections of rare and endangered oaks is a challenging task because they may frequently hybridize with other oak species. In other words, seedlings collected from a maternal tree in a botanic garden may be unintentional hybrid offspring that are not ideal for conservation purposes. This limits conservationists ability to use seeds produced in ex situ collections for restoration of threatened and endangered plant species unless the parentage is verified and the seedlings are not hybrids. 

We performed a study analyzing the parentage of acorns produced by maternal <i>Quercus muehlenbergii</i> individuals at the Morton Arboretum to identify if (1) any offspring individuals were hybrids and (2) what factors contribute to the parentage of offspring produced in living collections. We collected a total of 385 seeds, grew them to seedling stage, and sampled leaf tissue from the seedlings. Using 14 microsatellite loci, we performed parentage analysis using CERVUS software.

## Folder Descriptions

### Archive:
The Archive folder includes all files that are **not** integral to the analysis of this project. Inside it, there are three folders: Data_Files_Archive, UHA_Attempted_md_Analysis, and UHA_nomd_PolyPatEx_Analysis. The organization of the Archive folder is as follows:

  - Data_Files_Archive → This folder contains two duplicate data files that are not used in the final analysis.
  - UHA_Attempted_md_Analysis → This folder contains two folders: UHA_Cervus_Attempts and UHA_PolyPatEx_Attempts. The analysis included in these folders was done with files that had too much missing data so they are not relevant to the rest of our study. The purpose of the files was only for practicing the types of analysis with the data we had at the time. 
  - UHA_nomd_PolyPatEx_Analysis → This folder contains the fully cleaned data set analysis using PolyPatEx. For this project, using the CERVUS analysis gave us better, more accurate parentage assignment results, so we decided against using the PolyPatEx results. All the data and code files can be found here for PolyPatEx analysis.


### Data_Files:
The Data_Files folder includes two subfolers: Clean_Data_Files and Data_Cleaning_Code. The purpose of this folder is to include all the clean data files and to show how they were cleaned. These cleaned data files were used in the  paternity analysis to produce our final results.

  - Clean_Data_Files → This folder contains two subfolders, Clean_csv and Clean_xlsx. The files are the same, just in two different forms. These files are the results from the Data_Cleaning_Code folder.
  - Data_Cleaning_Code → The Data_Cleaning_Code folder contains an important script: UHA_datacleaning.R, written by Emily Schumacher. This R script cleans the data in order to make it useable for the analysis. The files directly involved with the data cleaning R script are "UHA_Final_Scores.arp", "UHA_Final_Scores.gen", and "UHA_Final_Scores_nomd_genalex.csv". The other csv files that were not directly involved with this script have minor adjustments made to them so they can be used in analysis as well. These files are "UHA_Final_Scores_nomd_df.csv" and "UHA_Final_Scores_GenAlEx.csv".

### Paternity_Analysis:
The Paternity_Analysis folder includes the cleaned data file used for analysis and code for the CERVUS analysis. The organization of the paternity analysis is as follows:

  - Results → This folder contains a R Markdown script written by Mikaely Evans that was used to create the final figures for this project. These figures were used in her symposium presentation, but are now outdated because this project is still in process. This folder also contains "UHA_parentage_results.csv" which is the csv file used to create these figures. It includes variables like Offspring ID, Mother ID, Candidate father ID, and Candidate father species.
  
  - UHA_nomd_CERVUS_Analysis → This folder contains the final no missing data parentage analysis input and output files from this study. CERVUS is a software made specifically for parentage analysis and it uses the three specific data sets mentioned in the CERVUS_Input_Files folder.
