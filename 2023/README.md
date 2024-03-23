# 2023
This folder contains scripts for discovering the governing equations of hydrodynamics from gridded data. This is the most recent version of SPIDER applied to the manuscript.

# Downloading Data
One design choice made during the rewrite of this paper was to not use the turbmat interface, and instead pull data from the cutout web service http://turbulence.idies.jhu.edu/cutout/ to access uninterpolated data directly. We did not want the interpolation process to introduce any uncertainty into our analysis. Downloading data in this way is not an automatic process, but the exact parameters to input into the web service are given in table 5 in appendix A. These h5 files need to be placed in a particular folder with particular names:

h5\_cutouts/velocity\_middle1.h5
h5\_cutouts/velocity\_edge1.h5
h5\_cutouts/pressure\_middle1.h5
h5\_cutouts/pressure\_edge1.h5


# Scripts
_analytic_test.m_ - computes an anyltic integral to show truncation error scaling for uniform and non-uniform grids.
_master_script.m_ - performs the weak-form integrals and performs sparse regression to obtain coefficients. Sections of this script should reproduce all figures of interest from the main body of the paper.
