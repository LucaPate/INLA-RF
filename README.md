# SPDE-Forest

## Scripts

### Strong non-linearity example

1. [runme_SPDE_strong_datasim.R](./Scripts/runme_SPDE_strong_datasim.R): Simulation of the data for the *strong non-linearity* scenario.
2. [runme_SPDE_strong_offset.R](./Scripts/runme_SPDE_strong_offset.R): Script implementing the INLA-RF algorithm to analyze the simulated data, correcting only the mean.
3. [runme_SPDE_strong_offsetUncertainty.R](./Scripts/runme_SPDE_strong_offsetUncertainty.R): Script implementing the INLA-RF algorithm to analyze the simulated data, correcting the mean and transfering the uncertainty between the steps of the algorithm.
4. [runme_SPDE_strong_lowRankNodeCorrection.R](./Scripts/runme_SPDE_strong_lowRankNodeCorrection.R): Script implementing the INLA-RF algorithm to analyze the simulated data, correcting only some pre-specified nodes of the spatio-temporal structure of the model.

### Other examples implementing the INLA-RF algorithm 

1. [runme_Data_analysis_MeuseData_INLARF_algorithm.R](./Scripts/runme_Data_analysis_MeuseData_INLARF_algorithm.R):
2. [runme_Data_analysis_Strong.R](./Scripts/runme_Data_analysis_Strong.R):
3. [runme_Jumps_StressPoints_INLARF_algorithm.R](./Scripts/runme_Jumps_StressPoints_INLARF_algorithm.R):

## Vignettes

1. [SPDE_RF_only_space.qmd](./Vignettes/SPDE_RF_only_space.qmd): The [html](./Vignettes/SPDE_RF_only_space.html) document. 
