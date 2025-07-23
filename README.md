# INLA-RF

## Scripts

### Simulation study 1 (spatio-temporal scenario)

1. [runme_spatiotemporal_simulation.R](./Scripts/runme_spatiotemporal_simulation.R): Script implementing the data simulation.
2. [runme_spatiotemporal_INLARF11.R](./Scripts/runme_spatiotemporal_INLARF11.R): Script implementing the INLA-RF1.1 algorithm to analyze the simulated data, correcting only the mean.
3. [runme_spatiotemporal_INLARF12.R](./Scripts/runme_spatiotemporal_INLARF12.R): Script implementing the INLA-RF1.2 algorithm to analyze the simulated data, correcting the mean and transfering the uncertainty between the steps of the algorithm.
4. [runme_spatiotemporal_INLARF11_and_INLARF12_CV.R](./Scripts/runme_spatiotemporal_INLARF11_and_INLARF12_CV.R): Script implementing the INLA-RF algorithms (1.1 and 1.2) to perform the cross-validation analysis of the simulated data, correcting the mean without and with uncertainty propagation.

### Simulation study 2 (purely temporal scenario)

1. [runme_temporal_INLARF2.R](./Scripts/runme_temporal_INLARF2.R): Script implementing the data simulation and the INLA-RF2 algorithm to analyze the simulated data, correcting only some pre-specified nodes of the spatio-temporal structure of the model.

### Extra

1. [diagonal_product.cpp](./Scripts/diagonal_product.cpp): Script defining a function to compute only the diagonal elements of a product of two square sparse matrtices in CSC format (for KLD computation).
