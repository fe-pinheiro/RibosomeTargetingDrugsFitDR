## RibosomeTargetingDrugsFitDR

Source code to fit dosage response curves of ribosome targeting drugs – Pinheiro et al., Nature Ecology & Evolution 2021. 

This folder contains custom MatLab scripts (v2016b) used to fit the metabolic model describing the dosage-response curves of *E. coli* exposed to the ribosome-targeting drug streptomycin as reported in Extended Data Fig. 1 of the paper. Optimization steps use fminsearch and lsqnonlin. We also use the parallel computing toolbox for speeding up computations whenever required. 


### Files

**1.	FitDR.m:** main script to fit the dosage-response curve of a given organism.

**2.	loadFunctionsGreulichModel:** solutions of the cubic equation that describes growth inhibition upon exposure to ribosome-targeting drugs for a given organism as described in Mol Syst Biol 2015 Mar;11(3):796. doi: 10.15252/msb.20145949

**3.	makeG:** prepare dosage-response for plotting

**4.	sseval3Pars.m:** function to fit the metabolic model with 3 parameters with fminsearch

**5.	sseval2ParsLs.m:** function to fit the constrained metabolic model (fixed Lambda\*) with fminsearch

**6.	sseval2ParsLsIC50.m:** function to fit the constrained metabolic model (fixed product of Lambda d\*) with lsqnonlin

**7.	getErrorBars.m:** 10-90 percentiles of basic parameters/observables

**8.	ensembleParsWT.mat and ensembleParsMM#.mat:** ensemble of parameters for the wild-type and membrane mutants 1-23 used to estimate bounds of response parameters/observables used in the paper. Files are provided for illustrative purposes and can be used in the example of getErrorBars.m

Further instructions, comments and examples are given in the files.

### Input data & running instructions

Dosage-response curves are fit to the average growth obtained from independent measurements of three replicates. Input data for the drug, growth values, and standard deviation of growth measurements are available from Supplementary Table 1. 

#### (i)	Run FitDR.m:

•	Choose of your favorite organism and define *drugData*, *growthData*, and *stdGrowthAux*. 

•	Run the first module to obtain optimal fit parameters and to plot the corresponding dosage-response curve. This module is used to generate all panels of Extended Data Fig. 1.

•	Run the second module to generate the input for the error analysis. Recommended use: save the variable *savePars* of the wild-type and of your favorite mutant to input them in getErrorBars.m. Error analysis in the manuscript was performed on an ensemble of parameters extracted from independent fits to 1000 synthetic dosage-response curves (accompanying .mat files).


#### (ii)	Run getErrorBars.m 

Input *saveParsWT* and *saveParsMut* generated by FitDR.m or use .mat files (see the commented example in getErrorBars.m). This script illustrates the propagation of error of scaled observables and gives the bounds of basic parameters/observables. This file can be easily adapted to compute bounds of more complex observables as desired. Analytical expressions are given in the paper. 

Please contact the authors if you have any specific question.
