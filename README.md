# Efficient Suzuki-Trotter decompositions

An implementation of a large variety of Suzuki-Trotter decomposition schemes (or splitting methods).

This repository contains the scripts and data required to reproduce the results presented in *"Optimised Trotter Decompositions for Classical and Quantum Computing"*, [arXiv:2211.02691 [quant-ph]](https://arxiv.org/abs/2211.02691), J. Phys. A: Math. Theor., DOI: [10.1088/1751-8121/acde7a](https://doi.org/10.1088/1751-8121/acde7a) as well as *"Simple Ways to improve Discrete Time Evolution"*, [arXiv:2309.xxxx [quant-ph]](https://arxiv.org/abs/2309.xxxxx), LATTICE 2023 proceedings.

For questions concerning the code contact [J.Ostmeyer@liverpool.ac.uk](mailto:J.Ostmeyer@liverpool.ac.uk).

## Derivations

The different schemes' theoretical efficiencies have been calculated in the Mathematica notebook `theo-efficiency.nb`. It also contains the complete code deriving the optimal 4th order decomposition schemes with real and complex coefficients.

Simply follow the examples if you want to calculate the efficiency of your own decomposition scheme.

The coefficients required for a factorised implementation of Taylor and Chebyshev series are calculated in `taylor_coefficients.nb`.

## Implementation

The different schemes are implemented in `C` with an `R` front-end. The code is located in the `simulation` directory.

The implementation can easily be augmented by another decomposition (see `ising_trotter.c`) and the Hamiltonian can be exchanged (see `ising_hamiltonian.c`) without any need to modify the Trotterizations.

### Compile
You might have to update the path to your `R/include` library in the `Makefile`. Find the path by executing `Sys.getenv("R_INCLUDE_DIR")` in your R environment.

After adjusting the Makefile according to your needs (possibly switching to another compiler), simply type `make` to compile all the `C` files.

### Run Code
The only high-level function `noisy_trace` calculates a trace (e.g. needed for the Frobenius norm) and can be executed by running `source("ising_dos.R")` in an `R` environment and then simply calling it.

The data used in the paper can be reproduced with the help of the `benchmark.R` script (e.g. run with `Rscript benchmark.R` in the console). Per default it is compiled without parallelisation and can take about an hour to run.

## Data and Plots

The data produced for the paper is located in `simulation/benchmark`.

The plots in the paper have been produced with `gnuplot`. The corresponding script `plot_tex.gp` can be found in `simulation/benchmark` as well.

## Stable Releases

`v1.0.0` initial publication with the preprint, compatible with arXiv versions v1, v2, v3.

`v1.1.0` added implementations of Yoshida's and Morales's methods, compatible with journal article and arXiv version v4.

`v1.2.0` added factorised polynomial implementations, compatible with LATTICE 2023 proceedings.
