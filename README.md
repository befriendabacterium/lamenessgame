README
=======

<p align="center">
<img src="https://github.com/befriendabacterium/lamenessgame/blob/main/game_screenshots.png" width="80%" height="80%">
</p>

Code for manuscript 'Exploring the potential of using simulation games for engaging with sheep farmers about lameness recognition'. The data are stored at https://osf.io/a6qu4/ and can be downloaded either manually (if you just want to look at the models/plots) or by running this pipeline (step 1 downloads the data from OSF). You can play the game here https://wheres-woolly.itch.io/lameness-game.

## Description

* author = Matt Lloyd Jones
* web = https://github.com/befriendabacterium/
* date = February 6th, 2023
* description = This repository contains a pipeline the reproduce the analyses in the preprint 'Exploring the potential of using simulation games for engaging with sheep farmers about lameness recognition' (https://www.biorxiv.org/content/10.1101/2022.10.26.513828v2).

## Pipeline

*N.B. Please use RStudio and load the 'lamenessgame_repo.Rproj' R Project file before running the scripts individually, in order. This will get you in the right directory straight away and allow you to dive straight in. You can run it in other ways and without RStudio, but it's easier this way.*

### 0. Acquire and load necessary packages

  * **script**: `0_acquirepackages.R` - Run this script to install and load the necessary packages for the pipeline. In brief, it goes through all the source code files, identifies packages called via the :: double colon operator, compares it to what you have installed, and installs and loads them if you don't.
  
  * **inputs**: None. 
  
    * **outputs**: Necessary packages that you didn't have before are downloaded to your R packages library.
  

### 1. Download data from OSF

  * **script**: `1_downloaddata.R`: Run this script to download the raw input data from OSF (https://osf.io/a6qu4/), plus some outputs that was produced manually (rather than via code).

  * **inputs**: None. 

  * **outputs**:

    * `inputs` folder: contains the raw data collected from participants in the study via MS_forms (`inputs/lamenessstudydata_010721.xlsx`), and an .RDS file with the desired order in which to rearrange the columns (`newcolnames.RDS`). 
    * `outputs` folder: folder that will contain the final data produced by the pipeline after processing and analysis of this data. Upon download, it only contains manually produced outputs from the study - Figure 1/game screenshots (`game_screenshots.csv`), the budget (`budget.csv`), and Supplementary Figure 1 of game screenshots relating to the game development process.
    
### 2. Process raw data

  * **script**: `2_wrangledata.R`: Run this script to wrangle the data into a format more suited to analysis and plotting. 
  
   * **inputs**: 
      * `inputs/lamenessstudydata_010721.xlsx`: Raw file of participant data outputted from MS forms.
      * `newcolnames.RDS`: R data file with desired order in which to rearrange columns
      
   * **outputs**:
      *  `processed_data/studydata_formatted.csv/RDS`: All participant performance data (i.e. participant feedback data form Likert questionnaire and open-form feedback has been removed) after wrangling and formatting via the pipeline. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      *  `processed_data/symptomslookedfor_df.csv/RDS`: A separate dataframe of the parsed (1 column per lameness sign/symptom) of the lameness signs/symptoms that each participant looked for data after wrangling and formatting via the pipeline. This data is already attached to `studydata_formatted`, but for the sake of plotting Figure 4 (`outputs/figures/symptoms.tiff`) it's easier to have it in a separate dataframe. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      *  `processed_data/symptomslookedfor_byfarmingexp_sum.csv/RDS`: A contingency table of the number of participants looking for each lameness sign, according to farming experience. This is used for the chi-square test of whether symptoms looked for differ according to farming experience, and the balloon plot, Supplementary Figure 3 (`outputs/figures/balloonplot.tiff`). The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      *  `processed_data/likertdata_formatted.csv/RDS` - Likert Questionnaire data after wrangling and formatting via the pipeline. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      
### 3. Analyse and plot data

  * **script**: `3_analysedata.R`: Run this script to download the data from OSF.

  * **inputs**:
    *  `processed_data/studydata_formatted.csv/RDS`: All participant performance data (i.e. participant feedback data form Likert questionnaire and open-form feedback has been removed) after wrangling and formatting via the pipeline. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
            *  `processed_data/symptomslookedfor_df.csv/RDS`: A separate dataframe of the parsed (1 column per lameness sign/symptom) of the lameness signs/symptoms that each participant looked for data after wrangling and formatting via the pipeline. This data is already attached to `studydata_formatted`, but for the sake of plotting Figure 4 (`outputs/figures/symptoms.tiff`) it's easier to have it in a separate dataframe. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
            *  `processed_data/symptomslookedfor_byfarmingexp_sum.csv/RDS`: A contingency table of the number of participants looking for each lameness sign, according to farming experience. This is used for the chi-square test of whether symptoms looked for differ according to farming experience, and the balloon plot, Supplementary Figure 3 (`outputs/figures/balloonplot.tiff`). The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      *  `processed_data/likertdata_formatted.csv/RDS` - Likert Questionnaire data after wrangling and formatting via the pipeline. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)

  * **outputs**:
     * `models/power_analysis.RDS`: Power analysis based on the sample size.
     * `models/NAME_model.RDS`: R data files of the linear models used in the analysis. The farming experience and user engagement models are built on the farming experience and time spent playing explanatory variables, respectively. The lameness signs/symptoms models are named with the convention 'symptom_symptomname_model'.
     * `models/NAME_adjp.RDS`: R data files of the adjusted p values linear models used in the analysis. Bonferroni correction was applied based on the number of previous models tested before that model was tested.
     * `models/symptomsVSfarmingexp_chisq.RDS`: Chi-squared test of whether lameness signs/symptoms participants looked for differed according to farming experience.
     * `figures/accuracyvsrecall.tiff`: Figure 2 in the manuscript, comparing accuracy and recall distributions for all participants (n=63)
     * `figures/farmingexperience.tiff`: Figure 3 in the manuscript, displaying relationships between recall scores and variables relating to participants' farming experience.
     * `figures/symptoms.tiff`: Figure 4 in the manuscript, displaying relationships between recall scores and variables relating to the 9 categories of lameness signs/symptoms looked for.
     * `figures/UE.tiff` - Figure 5 in the manuscript, displaying relationships between recall scores and variables relating to the user engagement with the game.
     * `figures/likertplot.tiff`: Figure 5 in the manuscript, a Likert plot of the Likert-Scale questionnaire data giving feedback on the game.
  
  
  
  

      
      
      
      
      

    
