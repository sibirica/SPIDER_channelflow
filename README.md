# SPIDER_channelflow
This is an open-source repository of the MATLAB code used in the paper Gurevich, Daniel R., Matthew R. Golden, Patrick A. K. Reinbold, and Roman O. Grigoriev, [*Learning fluid dynamics using sparse physics-informed discovery of empirical relations.*](https://arxiv.org/abs/2105.00048) arXiv preprint arXiv:2105.00048 (2021).

This code illustrates the sparse physics-informed discovery of empirical relations (SPIDER) algorithm, which is a general framework for robust identification of a complete model in terms of a system of partial differential equations and algebraic equations, comprising both bulk equations and boundary conditions. Input data representing a numerical simulation of a turbulent channel flow is taken from the Johns Hopkins University turbulence database: http://turbulence.pha.jhu.edu/Channel_Flow.aspx.

## Versions
Two versions of this code have been produced to discover governing equations from data. The earliest version of the manuscript uses the code from the 2021 directory, although it suffers from a misinterpretation of pressure data. All analysis has been redone and extended by independent code in the 2023 folder, which is the basis of the most recent manuscript. Both codebases are stored here for reproducibility of all variants of the manuscript.

There is a plan to continue development on a Python implementation of SPIDER, which can be found at https://github.com/sibirica/SPIDER_discrete.

