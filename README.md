# CADsimulation

# List of the functions

## CADsim.m

Function to perform the classical simulation of a noisy state ensemble, given a fixed set of unitaries for the simulation boxes.

## rho_d3.m & rho_d5.m

Ensembles {\rho_x}_{x=1}^{m}, with (d,m) = (3,5) and (d,m) = (5,3), used in the paper to evaluate the performance of the different optimization methods.

## Unitary randomization

It performs the classical simulation of a given ensemble, randomizing over the unitaries for the simulation.

## Unitary optimization

It performs the classical simulation of a given ensemble, optimizing over a fixed number of unitaries for the simulation.

## MUB rotations

When the ensemble has a specific structure, one can fix a starting set of simulation boxes and optimize over their unitary transformations.
In this case the starting set of simulation boxes is the set of unitaries associated to the MUBs in a d-dimensional Hilbert space, when d is odd.

# References:
G. Cobucci, A. Bernal, M. J. Renner, A. Tavakoli (2025). Operationally classical simulation of quantum states (https://arxiv.org/abs/2502.01440)

