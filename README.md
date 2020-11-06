# spectral_analysis_matrix_completion
Spectral analysis of time series data via subspace-based methods and  matrix completion based methods.

This repository contains classic subspace-based methods such as MUSIC [[S86]]

## List of algorithms

* `ANM` (_Atomic Norm Minimization_)  by author [[A20]](https://doi.org/10.1016/j.acha.2015.08.003).

# Installation
* Open MATLAB and run `setup_MatrixIRLS` or, alternatively, add subfolders manually to path. 

The main implementation can be found in the folder `algorithms/MatrixIRLS`, reference algorithms can be found in subfolders of `algorithms`. Some of the algorithms use methods for sparse access of factorized matrices implemented in C compiled as `.mex` files, and other auxiliary files contained in the `tools` folder.
## Examples and Usage
* `demo_MatrixIRLS`: This is a minimal example of how to use the main implementation  `MatrixIRLS.m`.

In the folder `examples/`, you find the following example scripts:
* `example_compare_MC_algos_badlycond`:
Compares several algorihms for low-rank matrix completion for the completion of badly conditioned matrix with condition number $\kappa=:\frac{\sigma_1(\textbf{X}_0)}{\sigma_r(\textbf{X}_0)} = 10^6$.  For this instance, it calculates 

## About this repository
##### Developer: 
* Thomas Weinberger (<t.weinberger@tum.de>)

If you have any problems or questions, please contact me for example via e-mail.

## References
 - [[S86]](https://ieeexplore.ieee.org/abstract/document/1143830) R. Schmidt, "Multiple emitter location and signal parameter estimation," in IEEE Transactions on Antennas and Propagation, vol. 34, no. 3, pp. 276-280, March 1986, doi: 10.1109/TAP.1986.1143830.
