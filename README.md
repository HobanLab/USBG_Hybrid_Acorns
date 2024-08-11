## Project Description
This GitHub repository is a combination of the analyses that were used to detect hybrids in the offspring of 9 _Quercus muehlenbergii_ trees sampled at the Morton Arboretum in Fall of 2022. These analyses were performed through 2022 - 2024 by a combination of collaborators: Ash Hamilton, Mikaely Evans, and Emily Schumacher. Our study was mostly concerned with quantifying the levels of hybridization within the botanic garden collections of Morton Arboretum white oaks to provide a model system for protecting oaks in living collections without producing hybrid offspring. 

We performed a study analyzing the parentage of acorns produced by maternal <i>Quercus muehlenbergii</i> individuals at the Morton Arboretum to identify if (1) any offspring individuals were hybrids and (2) what factors contribute to the parentage of offspring produced in living collections. We collected a total of 385 seeds, grew them to seedling stage, and sampled leaf tissue from the seedlings. Using 13 microsatellite loci, we performed parentage analysis using CERVUS software. Following this analysis, we produced multiple figures comparing the distance between candidate fathers, mothers, and hybrid production. 

## Project Workflow 
This project went through several stages, and so the workflow through the data files is described below: 

<b>Geographic analyses:</b> This was a preliminary stage of analyses that were performed to determine which botanic garden was a suitable model for our experimental question. These analyses were conducted using the RScripts stored in the "Geographic Analysis" file in the Analysis folder, with the data files used in the Data_Files/Geographic_Files pathway. The resulting CSV files are stored in the Results/Geographic_Analyses pathways.
- We examined 3 different botanic garden model systems: the UC Davis campus, Starhill Arboretum, and the Morton Arboretum. This stage of the process examined these arboretums for a few different traits:
    - At least 10 individuals of a specific oak species within breeding distance (350 m) that were producing acorns in the summer of 2022.
- Following the stage of this analysis we determined that the Morton Arboretum was the best suited to our study design.

<b>Parentage analyses</b>: Once the Morton Arboretum was selected as the study location, the geographic analyses at that site were used to generate sampling plans. _Quercus muehlenbergii_ was selected as the model species because it was the only oak species observed with 10 individuals over 10 years of age (reproductive age in oaks) within the study area that were producing at least 10 acorns. We sampled as many acorns as possible (with a goal of at least 10 acorns) and potting a minimum of 10 acorns. Acorns and leaf tissue for maternal individuals were sampled in Fall 2022, and then leaf tissue all white oak individuals within a 350 m radius from all maternal individuals were sampled in summer 2023 and summer 2024. DBH measurements were collected in summer 2024. Genetic data was collected from individuals in 2023 - 2024 which were used to perform parentage analysis in CERVUS. 
- CERVUS files: The RScripts stored in the Parentage_Analysis folder in the Analysis folders detail the data cleaning steps used to prep data files for parentage analysis, with data files stored in the Data_Files/CSV_Files folder or the Genotype_Files folder. The files generated in <u>01_data_cleaning_for_parentage.R"</u> were used to generate the input files for parentage analysis in CERVUS. The input files for CERVUS are stored in the <u>Data_Files/CERVUS_Files</u> pathway. There is a separate folder for <u>All_Loci</u> and <u>Red_Loci</u> because two forms of this analysis were performed - with _all_ 13 loci that were genotyped and then with a reduced set of loci - 9 - because these loci were identified to have a high frequency of null alleles (over 15%) - something that has been identified to bias parentage analyses. Within these parentage files there is a separate file for Input and Output files. Each run of parentage analysis in CERVUS requires an allele frequency file (contains allfreq within the title and .alf file), a simulation file for parentage analysis (includes sim in the name and is a .sim file), an offspring file (a CSV file that contains the word offspring in the title and is a CSV file), and a genotype file (titled genotype file with the word genotype_df in the title). These files are used to generate parentage results. The output files are the results of the CERVUS runs in the form of a text and CSV files.
- Parentage analysis result generation: The other RScripts in the <u>Parentage_Analysis</u> folder are to process the results of parentage analysis. The <u>02_data_cleaning_post_parentage.R</u> script is used to process the results of the CERVUS parentage runs - stored in the <i>par_sum</i> files - to produce figures and tables for the manuscript. This script includes also creating data files with the designation HCF. 
  
## Folder Descriptions

### Analysis:
- <b>Parentage_Analysis</b>: This folder contains all of the files used to perform parentage analysis with and without loci with a high amount of null alleles (all loci vs. red loci). The genotype files titled UHA_all_loci_genotype_df and UHA_red_loci_genotype_df are generated in the RScript "01_data_cleaning_for_parentage.R" in the RScripts folder. The offspring files are consistent between all and reduced loci praentage runs. The rest of the files in these folder are generated by the CERVUS parentage analysis software.
    - All_Loci
        - Input_Files
            - UHA_all_loci_genotype_df
            - all_loci_sim2.sim
            - all_loci_sim2.text
            - all_loci_af.alf    
        - Output_Files
            - all_loci_par_sum.csv
            - all_loci_par_sum.txt
    - Red_Loci
        - Input_Files
            - UHA_offspring2.csv
            - UHA_red_loci_genotype_df.csv
            - red_loci_allfreq.alf
            - red_loci_allfreq.txt
            - red_loci_sim.sim
            - red_loci_sim.text
        - Output_Files
            - red_loci_par_sum.csv
            - red_loci_par_sum.txt
           
- RScripts
    - 01_data_cleaning_for_parentage.R
    - 02_data_cleaning_post_parentage.R
    - 03_figures.R
    - 04_dist_analysis.R

### Archive:
The Archive folder includes files that were generated in the process of this 

  - Data_Files_Archive → This folder contains two duplicate data files that are not used in the final analysis.
  - UHA_Attempted_md_Analysis → This folder contains two folders: UHA_Cervus_Attempts and UHA_PolyPatEx_Attempts. The analysis included in these folders was done with files that had too much missing data so they are not relevant to the rest of our study. The purpose of the files was only for practicing the types of analysis with the data we had at the time. 
  - UHA_nomd_PolyPatEx_Analysis → This folder contains the fully cleaned data set analysis using PolyPatEx. For this project, using the CERVUS analysis gave us better, more accurate parentage assignment results, so we decided against using the PolyPatEx results. All the data and code files can be found here for PolyPatEx analysis.


### Data_Files:
The Data_Files folder includes two subfolers: Clean_Data_Files and Data_Cleaning_Code. The purpose of this folder is to include all the clean data files and to show how they were cleaned. These cleaned data files were used in the  paternity analysis to produce our final results.

  - Clean_Data_Files → This folder contains two subfolders, Clean_csv and Clean_xlsx. The files are the same, just in two different forms. These files are the results from the Data_Cleaning_Code folder.
  - Data_Cleaning_Code → The Data_Cleaning_Code folder contains an important script: UHA_datacleaning.R, written by Emily Schumacher. This R script cleans the data in order to make it useable for the analysis. The files directly involved with the data cleaning R script are "UHA_Final_Scores.arp", "UHA_Final_Scores.gen", and "UHA_Final_Scores_nomd_genalex.csv". The other csv files that were not directly involved with this script have minor adjustments made to them so they can be used in analysis as well. These files are "UHA_Final_Scores_nomd_df.csv" and "UHA_Final_Scores_GenAlEx.csv".

### Results:
