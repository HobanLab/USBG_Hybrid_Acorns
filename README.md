## Project Description
This GitHub repository is a combination of the analyses that were used to detect hybrids in the offspring of 9 _Quercus muehlenbergii_ trees sampled at the Morton Arboretum in Fall of 2022. These analyses were performed through 2022 - 2024 by a combination of collaborators: Ash Hamilton, Mikaely Evans, and Emily Schumacher. Our study was mostly concerned with quantifying the levels of hybridization within the botanic garden collections of Morton Arboretum white oaks to provide a model system for protecting oaks in living collections without producing hybrid offspring. 

We performed a study analyzing the parentage of acorns produced by maternal <i>Quercus muehlenbergii</i> individuals at the Morton Arboretum to identify if (1) any offspring individuals were hybrids and (2) what factors contribute to the parentage of offspring produced in living collections. We collected a total of 385 seeds, grew them to seedling stage, and sampled leaf tissue from the seedlings. Using 13 microsatellite loci, we performed parentage analysis using CERVUS software. Following this analysis, we produced multiple figures comparing the distance between candidate fathers, mothers, and hybrid production. 

## Project Workflow 
This project went through several stages, and so the workflow through the data files is described below: 

<b>Geographic analyses:</b> This was a preliminary stage of analyses that were performed to determine which botanic garden was a suitable model for our experimental question. These analyses were conducted using the RScripts stored in the "Geographic Analysis" file in the Analysis folder, with the data files used in the <b>Data_Files/Geographic_Files</b> pathway. The resulting CSV files are stored in the <b>Results/Geographic_Analyses</b> pathways.
- We examined 3 different botanic garden model systems: the UC Davis campus, Starhill Arboretum, and the Morton Arboretum. This stage of the process examined these arboretums for a few different traits:
    - At least 10 individuals of a specific oak species within breeding distance (350 m) that were producing acorns in the summer of 2022.
- Following the stage of this analysis we determined that the Morton Arboretum was the best suited to our study design.

<b>Parentage analyses</b>: Following our decision that the Morton Arboretum was the best site for this analysis, we sampled _Quercus muehlenbergii_ individuals for acorns and leaf tissue. The acorns were grown into seedlings and then leaf tissue was sampled for parentage analyses, and these individuals, as well as all maternal individuals and any white oak within a 350 m radius of each maternal individual was sampled for parentage analysis. These individuals were genotyped using 13 microsatellite loci and analyzed in CERVUS parentage software. The description of the data processing to prep files for parentage analysis and the processing of the outputs is described below. All of the files used for running parentage analysis and storing parentage analysis data files are in stored in the <b>Analysis/Parentage_Analysis</b> folder, with the analyses being run from the <b>RScripts</b> folder, and all of the data files stored in the <b>CERVUS_Files</b> folder. The results of parentage analyses are stored in the <b>Results/Parentage_Analysis</b> folder. 
- Preparing files for CERVUS: Genotype files generated from Geneious were processed in the <b>Analysis/Parentage_Analysis/RScripts /01_data_cleaning_for_parentage.R</b> file and then were used to generate the input files for parentage analysis in CERVUS. The input files used in this analysis were a genepop file (.genepop) and CSV file generated from Geneious files. These files were used to generate "clean" genotype files, which removed any individuals with 25% missing genotypes. We also tested loci for linkage disequilibrium and null allele frequency. Following null allele analysis, we identified 4 loci with high frequencies of null alleles (>15%) and so we created data files with and without these loci in this data file. "all_loci" files refer to data files with all 13 loci run on individuals in these data, whereas "red_loci" files refer to genotype files that have 4 loci removed because they were identifed to have high frequencies of null alleles (>15%). The results of this script are stored in the <b>Analysis/CERVUS_Files</b> pathway.  
- CERVUS_Files: All files used to run parentage analysis in CERVUS are stored in the <b>Analysis/CERVUS_Files</b> folder, which has separate folders for "all_loci" and "red_loci" data files, as these parentage runs were done separately. Each scenario of the analysis has an "Input_Files" and "Output_Files" folder. The "Input_Files" folder store all of the files needed to run parenatge analysis in CERVUS: a cleaned, genotype_df CSV file, an offspring file, an allele frequency file (.alf) generated in CERVUS, and a simulation file (.sim) generated in CERVUS. The genotype data file is identical to the cleaned score genotype file generated by the <b>01_data_cleaning_for_parentage.R</b> RScript, with all loci or reduced loci. The output folder contains the resulting parentage assignment CSV files with the suffix "par_sum". 
- Parentage analysis result generation: The other RScripts in the <u>Parentage_Analysis</u> folder are to process the results of parentage analysis. The <u>02_data_cleaning_post_parentage.R</u> script is used to process the results of the CERVUS parentage runs - stored in the <i>par_sum</i> files - to produce figures and tables for the manuscript. This script includes also creating data files with the designation HCF. HCF stands for "high confidence father" which are candidate father assignments that were made with pairwise and trio LOD scores > 0. This resulted in four <u>par_sum</u> data files - overall four scenarios to determine the impact of null alleles and confidence of parentage assignments on the final figures. This RScript was initially written by Mikaely Evans and Ash Hamilton but then udpated to loop over all data file scenarios by Emily Schumacher.
    - Final figures for the manuscript were generated in the <u>03_figures.R</u> RScript. This code was generated by Mikaely Evans.  
  
## Folder Descriptions

### Analysis:
This folder is divided into 3 separate analysis sections that cover a different set of analyses. 
- Geographic_Analysis
    - 01_geographic_analysis_prep.R
        - Description: This R Script was used to visualize oak species in botanic gardens Starhill Arboretum, UC Davis Campus, and the Morton Arboreutm. These data files were used to generate color coded maps stored in the project guide for this project. 
    - 02_garden_summary_dfs.R
        - Description: This R Script was used to generate overview data frames of the oak species in Starhill Arboretum, UC Davis Campus, and the Morton Arboreutm - oak individuals above 10 years of age to visualize candidate sampling areas.    
    - 03_TMA_all_trees.R
        - Description: This script was used to visualize candidate sites for acorn sampling once the Morton Arboretum was decided to be the best site for this project.
- Parentage_analysis
    - 01_data_cleaning_for_parentage.R
    - 02_data_cleaning_post_parentage.R
    - 03_figures.R
- STRUCTURE
    
  
### Archive:

### Data_Files:

### Results:
