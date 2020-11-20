# spectralAna_matComp
Spectral analysis of time series data via subspace-based methods and matrix completion.

This repository contains classic subspace-based methods such as MUSIC [[S86]], ESPRIT [[R+89]] and the Matrix Pencil Method [[H+90]] for harmonic retrieval from complete and equispaced samples as well as modern techniques from matrix completion such as Atomic Norm Minimization, Nuclear Norm Minimization and Iteratively Reweighted Least Squares, which can be used as a pre-processing step when the available data is incomplete. The performance of the different algorithms can be compared by running the included scripts, which give empirical results based on several Monte-Carlo simulations in which different parameters of the estimation problem are varied, respectively.

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
 - [[R+89]] (https://ieeexplore.ieee.org/document/32276) R. Roy and T. Kailath, "ESPRIT-estimation of signal parameters via rotational invariance techniques," in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 37, no. 7, pp. 984-995, July 1989, doi: 10.1109/29.32276.
 - [[H+90]] (https://ieeexplore.ieee.org/document/56027) Y. Hua and T. K. Sarkar, "Matrix pencil method for estimating parameters of exponentially damped/undamped sinusoids in noise," in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 38, no. 5, pp. 814-824, May 1990, doi: 10.1109/29.56027.
