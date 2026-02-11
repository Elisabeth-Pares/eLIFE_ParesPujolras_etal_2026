# eLIFE_ParesPujolras_etal_2026

This repository contains analysis scripts for the following paper:

[Parés-Pujolràs, E., Geuzebroek, A. C., O’Connell, R. G., & Kelly, S. P. (2025). Perceptual glimpses are locally accumulated and globally maintained at distinct processing levels. eLife, 14.
](https://elifesciences.org/reviewed-preprints/107980)

The scripts included in the repository and the data on OSF (COMING SOON) enable the replication of the data analysis and figures reported in the paper. All are shared under a CC-BY 4.0 license.

Code shared here was developed using Matlab R2021a, the [EEGab toollbox](https://sccn.ucsd.edu/eeglab/) (v.2021.0), [FieldTrip toolbox](https://www.fieldtriptoolbox.org/), [CSD toolbox](https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/), and the [TBT toolbox](https://github.com/mattansb/TBT).

## Task
Contains the scripts used to run both experiments in the paper. The experimental task is coded in Matlab, using Psychtoolbox-3 (http://psychtoolbox.org/). 
The folder includes:
- EP1_GapsTask.m: master script to run the experimental tasks in the paper. 
  - functions/: includes all required functions to run the task.
  - functions/ShadlenDotsx/: basic random dot motion task, retrieved from (https://shadlenlab.columbia.edu/resources/VCRDM.html) with minor edits.

## Analysis
Contains the scripts used to perform all analyses reported in the paper. Main contents:

- GT_Analysis_Paper.R:
  - Preprocesses behavioural data and produces the aggregate files loaded in the Matlab scripts plotting behavioural results ("EPX_procBehData.csv").
  - Statistical tests reported in Eq. 1 to 4.
  
- GT_EEGAnalysis.m: master script for EEG analysis. All required functions are included in the /EEG preprocessing/ folder. 
  - /modelling/:includes all required functions to fit the models in the paper, and to produce the behavioural/EEG simulations.
  - /figures/: reproduces Figures 1-5 in the main manuscript. Supplementary figures can be produced through this code too, by changing the appropriate parameters in each function and as commented in the code.
