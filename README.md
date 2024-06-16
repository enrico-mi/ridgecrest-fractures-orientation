### Code for _Coseismic damage of the 2019 Ridgecrest earthquake consistent with Mohr-Coulomb failure_

This repository contains the code developed for the analyses of the fractures' orientation published in the manuscript _Coseismic damage of the 2019 Ridgecrest earthquake consistent with Mohr-Coulomb failure_.

#### Requirements

The code was developed in the following environment:
* MATLAB Version 9.14 (R2023a)
* Mapping Toolbox Version 5.5 (R2023a)
* Parallel Computing Toolbox Version 7.8 (R2023a)
<!-- * Signal Processing Toolbox Version 9.2 (R2023a) -->
<!-- * Statistics and Machine Learning Toolbox               Version 12.5        (R2023a) -->

The Parallel Computing Toolbox is required only for the `parfor` loops in the `batch` files; the same files can easily be run serially by replacing those loops with classical `for` loops.

#### Instructions

###### Analyses

The analysis is organized in three stages:
1. compute the stresses at the fractures' centers
2. estimate errors between observed and predicted orientations
3. find best fits and estimate uncertainty

Each stage can be run with its corresponding "batch" file:
1. `batch_stresses.m`
2. `batch_errors.m` (`batch_errors_distance.m` for distance-dependent analysis)
3. `batch_minima.m`

These files are set up to automatically perform all the analyses at each stage, and _should_ work straight out of the box. To achieve this, all input/output is hardcoded to some degree, and the repository comes with empty folders (`output` and `output/angles_analysis`) where the outputs are stored.

Note that `batch_stresses.m` calls `run_calc_stresses_NW.m`, which computes the stresses at the fractures' centers by adding regional lithostatic stresses to previously computed coseismic stresses. It requires input files for the coseismic stresses, provided in the folder `cs-output`. If you want to use your own stress field, follow the structure of those files. It's a tabular format where each row corresponds to the stress computed at the center point of a given fracture, ordered as in the Ponti et al (2019) shape file. You'll have to replace replace your own filenames in `run_calc_stresses_NW.m`.

###### Plots

`SRw_plots_theta_H__NW` plots SRw as a function of the regional stress orientation in the absence of coseismic stresses (k=0, as Fig 2a in the manuscript).

`SRw_plots_magnfactor_NW` plots SRw as a function of parameter _k_ for a prescribed regional stress orientation (N14E, as Fig 2b in the manuscript).
