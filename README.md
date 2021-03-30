# SPIDER_channelflow
This is an open-source repository of the MATLAB code used in the paper Gurevich, Daniel R., Patrick A. K. Reinbold, and Roman O. Grigoriev, *Physics-informed learning of complete models from highly turbulent data.* arXiv preprint arXiv:xxxx.xxxxx (2021) [arXiv link here]

This code illustrates the sparse physics-informed discovery of empirical relations (SPIDER) algorithm, which is a general framework for robust identification of a complete model in terms of a system of partial differential equations and algebraic equations, comprising both bulk equations and boundary conditions. Input data representing a numerical simulation of a turbulent channel flow is taken from the Johns Hopkins University turbulence database: http://turbulence.pha.jhu.edu/Channel_Flow.aspx.

## Running the code
To get started, install Turbmat from https://github.com/idies/turbmat in order to be able to query the JHU turbulence database. Then, you may execute the main script *Run_SPIDER_NS.mat* to begin downloading the data and run the sparse regression from the Navier-Stokes equation (this may take some hours the first time before the data has been permanently saved locally). The purpose of each .mat script is explained below:

**Run_SPIDER_NS.mat** Main script for sparse regression, specifying various options in the algorithm. In particular, changing *mode* between "NS", "div", and "BC" allows for identification of the rank-1 and rank-0 equations and the boundary conditions, respectively.

**SPIDER_NS.mat** Routine for loading the data and integrating the library terms over several domains with respect to the specified weight functions.

**SparseReg.mat** Sparse regression routine identifying best single-term and multi-term models.

**glossaryNS.mat** Lists the libraries of terms used in the sparse regression for each equation -- some extra terms that were not used in the paper are commented out.

**defineLibrary.mat** Defines the library terms after integration by parts in order to reduce discretization error in *SPIDER_NS.mat*.

**diff_dim.mat** Helper function for computing central finite differences of boundary terms.

**pullPressure.mat, pullVelocity.mat** Helper functions for loading velocity and pressure data from the JHU database.

**LegendrePoly.mat, weight_3d.mat, etc.** Helper functions for computing Legendre polynomials and building weight functions over an integration domain.

Finally, *y.txt* contains the grid of y-coordinates of the turbulence data.
