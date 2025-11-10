# ADHD_deviation
Data and codes for our paper "Developmental deviations of structural connectivity in youths with ADHD predict symptom and treatment outcomes".  
We investigated how structural connectivity alterations unfold during youth and their clinical significance across two large cohorts: the U.S. [Adolescent Brain Cognitive Development study](https://abcdstudy.org/) (6,687 typically developing and 1,114 ADHD scans) and an independent Chinese cohort (292 typically developing and 468 ADHD participants). The Chinese cohort consist of two studies, EFNY (Executive Function and Neurodevelopment in Youth) and PKU6 (Peking University Sixth Hospital). Data of these two datasets are not yet available, as collection is still ongoing.  

## Software and system requirements
### Diffusion & structural MRI preprocessing
* FreeSurfer v7.1.1 (https://surfer.nmr.mgh.harvard.edu/)
* QSIPrep 0.16.0 (https://qsiprep.readthedocs.io/)
* OS: Linux

### Postprocessing
* Connectome Workbench v2.0.1 (https://www.humanconnectome.org/software/connectome-workbench)
* R v4.2.2 (https://www.r-project.org)
* MATLAB R2020a (https://www.mathworks.com/)
* OS: Windows / Linux

The system requirements and installation guide for each software can be found on its respective website.

## demography
This folder contains demographic tables for samples from the three studies. 

## Data
`centiles_Yeo17_*.rds`
Structure: a list of length 120. Each element contains two 3 × 1000 matrices:
1. Males: rows correspond to P2.5, P50, and P97.5; columns index age points.
2. Females: rows correspond to P2.5, P50, and P97.5; columns index age points.
Age grid: the 1,000 columns represent age points sampled uniformly across the dataset’s age range.

`MultiMediation_SCDeviation_*_Age_CV75_Yeo17_ADHDall_controlSex.rds`
Results of the mediation analyses.

`plotdatasum.df_Yeo17_sumSCdeviation_CV75_*.rds`
Data for visualizing developmental trajectories of deviations in structural connectivity strength.

`Yeo17_SArank.csv`
Median sensorimotor–association (S–A) rank of vertices within each Yeo network.

`schaefer376_index_Yeo17.csv`
Indices and network labels for the Schaefer-400 atlas with limbic regions removed (376 regions).

`SCdeviation.PCA.rds`
Principal component analysis (PCA) results for the Chinese cohort.

## ADHD_script
`S0_NetworkPlot`
Scripts for:
- plotting the S–A rank distributions across Yeo networks,
- generating the S–A connectional-axis matrix, and
- creating Yeo-17 network CIFTI files ordered along the S–A axis.

`S1_dataclean_merge`  
Scripts for screening participants and organizing structural-connectivity strengths for 120 connections into `N_obs × N_conn` data frames.  

`S2_normativemodeling`  
Scripts for constructing normative models and computing individual deviations per connection.  
[1] `S1_selectparameters_*.qmd`: Code to select the optimal distribution family and spline parameters for all connections.  
[2] `S2_bootstrap_*.R`: Defines a function that performs one iteration of the bootstrap analysis. Called by `S2_bootstrap_*_exe.R`.  
[3] `S3_Evaluate_GAMLSS_*.R`: Code for evaluating the normative models.  
[4] `S4_constructNM_forDeviation_*.R`: Constructs the normative models and computes individual deviations.  
[5] `S5_plotNM_eachGroup_*.qmd`: Visualizes developmental trajectories derived from the normative models and their key features; also plots mean deviations for both groups.  
[6] `S6_demodescription.R`: Generates age distributions and tables describing demographic information.  

`S3_directcompare`  
Scripts for comparing developmental trajectories between ADHD and TD groups and testing age effects on deviations.    
[1] `S1_compareSCDeviation_ADHDall_VS_TD_ABCD.qmd`: Performs an age × diagnosis interaction analysis on SC deviations.  
[2] `S2_DeviationChangeByAge_*.Rmd`: Tests age effects on SC deviations within each diagnosis group.  

`S4_association_ADHDsymp_SCstrength`  
Scripts for examining symptom relevance and treatment associations with deviations in structural connectivity strength.  
[1] `S1_ADHDsymp_SCstrengthDeviation_ABCD.Rmd`: Tests associations between SC deviations and ADHD symptoms in the ABCD dataset.  
[2] `S1_ADHDsymp_SCstrengthDeviation_Chinese.Rmd`: Tests associations between SC deviations and ADHD symptoms, and examines baseline deviations versus treatment responses in the Chinese cohort.  
[3] `S2_ADHDsymp_SCDeviation_Med_ABCD.qmd`: Performs mediation analyses of age, deviations, and symptoms.  
[4] `S3_delta_deviation_delta_symptom_ABCD.qmd`: Tests within-subject associations between rates of change in deviations and rates of change in symptoms.  
[5] `S4_ADHDsymp_medication_PKU6.qmd`: Examines changes in deviations following medication treatment.  

## functions & gamfunction
These two folders contain the R functions required for conducting the analyses in this study.  

