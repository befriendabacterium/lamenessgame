README
=======

<p align="center">
<img src="https://github.com/befriendabacterium/lamenessgame/blob/main/game_screenshots.png" width="80%" height="80%">
</p>

This repository contains the code associated with the manuscript 'Exploring the potential of using simulation games for engaging with sheep farmers about lameness recognition'. The data are stored at https://osf.io/a6qu4/ and downloaded, processed and analysed in this pipeline. To run the pipeline and reproduce the analysis and manuscript itself, please download (green 'Code' button, top right of this page) or clone (https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) this repository  to your local computer, open the R project file ('lamenessgame_repo.Rproj'), and then run each of the R scripts in 'src' in RStudio individually and in order (see 'Pipeline' below). 

If you want to see the outputs of the code without running it, then just head over to https://doi.org/10.5281/zenodo.7605244, where all code, input and output data are archived. If you're reading this having downloaded this Zenodo repository but want to run the code anyway, then you can skip straight to step 2 ('2. Process raw data') in the pipeline.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7605244.svg)](https://doi.org/10.5281/zenodo.7605244)

If you'd just like to play the game, then head over to https://wheres-woolly.itch.io/lameness-game to play it online, or to https://doi.org/10.5281/zenodo.7612059 to play it offline.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7612059.svg)](https://doi.org/10.5281/zenodo.7612059)


## Description

* author = Matt Lloyd Jones
* web = https://github.com/befriendabacterium/lamenessgame
* date = February 14th, 2023
* description = This repository contains a pipeline to reproduce the analyses in the preprint 'Exploring the potential of using simulation games for engaging with sheep farmers about lameness recognition' (https://www.biorxiv.org/content/10.1101/2022.10.26.513828v2).

## Pipeline

Before beginning anything, please open the 'lamenessgame_repo.Rproj' R Project file in RStudio before running the scripts individually, in the order described in the pipeline. This will get you in the right directory straight away and allow you to dive straight in. You can run it in other ways and without RStudio, but it's easiest this way.

### 0. Acquire and load necessary packages

  * **script**: `src/0_acquirepackages.R` - Run this script to install and load the necessary packages for the pipeline. In brief, it goes through all the source code files, identifies packages called via the :: double colon operator, compares it to what you have installed and installs them if you don't. Finally, it loads all the required packages via library() - **so you need to run this script to make the pipeline work even if you have all the required packages (because they may not be loaded)**
  
  * **inputs**: None. 
  
  * **outputs**: Necessary packages that you didn't have before are downloaded to your R packages library.

### 1. Download data from OSF

  * **script**: `src/1_downloaddata.R`: Run this script to download the raw input data from OSF (https://osf.io/a6qu4/), plus some outputs that was produced manually (rather than via code).

  * **inputs**: None. 

  * **outputs**:

    * `inputs` folder: contains the raw data collected from participants in the study via MS_forms (`inputs/lamenessstudydata_010721.xlsx`), and an .RDS file with the desired order in which to rearrange the columns (`newcolnames.RDS`). 
    * `outputs` folder: folder that will contain the final data produced by the pipeline after processing and analysis of this data. Upon download, it only contains manually produced outputs from the study - Figure 1/game screenshots (`game_screenshots.csv`), the budget (`budget.csv`), and Supplementary Figure 1 of game screenshots relating to the game development process.
       
### 2. Process raw data

 *N.B. If you've downloaded the Zenodo repository rather than downloaded/cloned the Github repo, you can start from here as you'll already have the input data (and the output data, actually).*

  * **script**: `src/2_wrangledata.R`: Run this script to wrangle the data into a format more suited to analysis and plotting. 
  
   * **inputs**: 
      * `inputs/lamenessstudydata_010721.xlsx`: Raw file of participant data outputted from MS forms.
      * `newcolnames.RDS`: R data file with desired order in which to rearrange columns
      
   * **outputs**:
      *  `processed_data/studydata_formatted.csv/RDS`: All participant performance data (i.e. participant feedback data form Likert questionnaire and open-form feedback has been removed) after wrangling and formatting via the pipeline. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      *  `processed_data/symptomslookedfor_df.csv/RDS`: A separate dataframe of the parsed (1 column per lameness sign/symptom) of the lameness signs/symptoms that each participant looked for data after wrangling and formatting via the pipeline. This data is already attached to `studydata_formatted`, but for the sake of plotting Figure 4 (`outputs/figures/symptoms.tiff`) it's easier to have it in a separate dataframe. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      *  `processed_data/symptomslookedfor_byfarmingexp_sum.csv/RDS`: A contingency table of the number of participants looking for each lameness sign, according to farming experience. This is used for the chi-square test of whether symptoms looked for differ according to farming experience, and the balloon plot, Supplementary Figure 3 (`outputs/figures/balloonplot.tiff`). The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      *  `processed_data/likertdata_formatted.csv/RDS` - Likert Questionnaire data after wrangling and formatting via the pipeline. The .csv format is more human-readable, whilst the .RDS file is more computer-readable and the input used in the pipeline, because it preserves all changes made in R needed to reproduce the outputs exactly (e.g. ordering of factor levels in plots)
      
### 3. Analyse and plot data

  * **script**: `src/3_analysedata.R`:  Run this script to run the quantitative analyses and (re)produce the R-produced figures used in the manuscript.

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
  
  
 ### 4. Write manuscript via R Markdown
 
 *N.B. Running R Markdown scripts works a little differently to running R scripts - open them in RStudio and then click the 'Knit' button (ball of twine with needle through it) instead of running R code as you usually would. You'll need 'rmarkdown' and a few other packages to get this to work - I've tried to include these in '0_acquirepackages.R' but if RStudio prompts you for more, then just click 'Yes' to install the packages. Another thing to note is that the outputted Word document may prompt you when opened with 'This document contains references to other files. Do you want to update the fields in this document - Yes/No' - just click 'Yes' and the Word doc should open fine.*
   
   * **script**: `lamenessgame_MS.Rmd`: Run this R Markdown script to reproduce the submitted manuscript. Note this sits outside of the 'src' folder, in the main directory.
   * **inputs**: 
      * `outputs`: Contents of the 'outputs' folder and its sub-folders, which contain all processed data, models and figures used in the paper.
      * `article_template.docx`: A template Word article which contains the necessary Styles to make the manuscript look how we want.
      * `lamenessgame_refs.bib`: Bibliography file containing the references used in the manuscript.
      * `harvard-cite-them-right_12thed_no-et-al.csl` - Citation style file to produce in-text citations and a bibliography in Cite Them Right 12th edition's Harvard style (no 'et al' used in the bibliography). Obtained from Zotero Style repository (https://www.zotero.org/styles?q=id%3Aharvard-cite-them-right-no-et-al) 
   * **outputs**: 
      * `lamenessgame_MS.docx` - A Word file of the submitted manuscript.
  

      
      
      
      
      

    
