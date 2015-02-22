# euler2d_cudaFortran
2nd order Godunov solver for 2d Euler equations written in CUDA Fortran

# Use CUDA Fortran (PGI compiler)

This code solves the 2D Euler equations in a regular cartesian mesh
using a 2nd order godunov-based finite volume scheme.

Parameters:
	edit file test.nml

There are 2 differents variants of the numerical scheme which are different in term of allocated memory. Use parameter implementationVersion = 0 or 1 to change.

Initial condition: a discontinuity along the domain diagonal

Example of use:
./euler2d_gpu ./test.nml

Output: 
	VTK ascii file using VTK Image Data format
Visualization: 
	paraview --data=euler2d_..vti


Pierre Kestener - October 12, 2013