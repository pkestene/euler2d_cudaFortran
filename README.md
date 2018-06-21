# euler2d_cudaFortran
2nd order Godunov solver for 2d Euler equations written in CUDA Fortran

WARNING: this code is deprecated. You may be interested to a new version using C++/kokkos for performance portability:
https://github.com/pkestene/euler2d_kokkos

## Short description

This code solves the 2D Euler equations in a regular cartesian mesh
using a 2nd order godunov-based finite volume scheme.

## PGI compiler for CUDA Fortran
It is written in CUDA Fortran and designed as a simple example of use of NVIDIA GPU's for CFD applications. You need to have the PGI compiler to build this application.

## Parameter file
Parameters:
	edit file test.nml

There are 2 differents variants of the numerical scheme which are different in term of allocated memory. Use parameter implementationVersion = 0 or 1 to change.

Initial condition: a discontinuity along the domain diagonal


## Example of use

./euler2d_gpu ./test.nml

Output: 
	VTK ascii file using VTK Image Data format

Visualization: 
	paraview --data=euler2d_..vti


Pierre Kestener - October 12, 2013
