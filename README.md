# INLA-RF

## Scripts

### Strong non-linearity example

1. [runme_spatiotemporal_simulation.R](./Scripts/runme_spatiotemporal_simulation.R): Simulation of the data for the *strong non-linearity* scenario.
2. [runme_spatiotemporal_INLARF11.R](./Scripts/runme_spatiotemporal_INLARF11.R): Script implementing the INLA-RF algorithm to analyze the simulated data, correcting only the mean.
3. [runme_spatiotemporal_INLARF12.R](./Scripts/runme_spatiotemporal_INLARF12.R): Script implementing the INLA-RF algorithm to analyze the simulated data, correcting the mean and transfering the uncertainty between the steps of the algorithm.
4. [runme_spatiotemporal_INLARF11_and_INLARF12_CV.R](./Scripts/runme_spatiotemporal_INLARF11_and_INLARF12_CV.R): Script implementing the INLA-RF algorithm to perform the cross-validation analysis of the simulated data, correcting the mean with and without uncertainty propagation..
5. [runme_temporal_INLARF2.R](./Scripts/runme_temporal_INLARF2.R): Script implementing the INLA-RF algorithm to analyze the simulated data, correcting only some pre-specified nodes of the spatio-temporal structure of the model. 

### Other examples implementing the INLA-RF algorithm 

1. [runme_Data_analysis_MeuseData_INLARF_algorithm.R](./Scripts/runme_Data_analysis_MeuseData_INLARF_algorithm.R):
2. [runme_Data_analysis_Strong.R](./Scripts/runme_Data_analysis_Strong.R):
3. [runme_Jumps_StressPoints_INLARF_algorithm.R](./Scripts/runme_Jumps_StressPoints_INLARF_algorithm.R):

## Vignettes

1. [SPDE_RF_only_space.qmd](./Vignettes/SPDE_RF_only_space.qmd): The [html](./Vignettes/SPDE_RF_only_space.html) document. 
