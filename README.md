# CADsimulation

# List of the functions

## CADsim

Function to perform the classical simulation of a noisy state ensemble, given a fixed set of unitaries for the simulation boxes.

## Unitary randomization

It performs the classical simulation of a given ensemble, randomizing over the unitaries for the simulation.

## Unitary optimization

It performs the classical simulation of a given ensemble, optimizing over a fixed number of unitaries for the simulation.

## MUB rotations

When the ensemble has a specific structure, one can fix a starting set of simulation boxes and optimize over their unitary transformations.
In this case the starting set of simulation boxes is the set of unitaries associated to the MUBs in a d-dimensional Hilbert space, when d is odd.
