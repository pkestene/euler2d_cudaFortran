# euler2d_cudaFortran

2nd order Godunov solver for 2d Euler equations written in CUDA Fortran

Note: If your more interested into c++ for GPUs, have a look at version using C++/kokkos for performance portability:
https://github.com/pkestene/euler2d_kokkos

## Short description

This code solves the 2D Euler equations in a regular cartesian mesh
using a 2nd order godunov-based finite volume scheme.

There are actually two version, one is using the so called [cuda fortran](https://developer.nvidia.com/cuda-fortran) programming model, and the other is designed with the newer [fortran standard parallelism](https://developer.nvidia.com/blog/accelerating-fortran-do-concurrent-with-gpus-and-the-nvidia-hpc-sdk/) (stdpar).

## Nvhpc compiler for Cuda Fortran and stdpar

You need to have installed [nvhpc] to build this application.

Build for a given Cuda architecture, using a given Cuda runtime version (assuming you install nvhpc)

```shell
# for the cuda fortran version
cd cuf
make CUDA_ARCH=cc80 CUDA_VERSION=11.6
# for the stdpar version
cd stdpar
make CUDA_ARCH=cc80 CUDA_VERSION=11.6
```

## Parameter file

Parameters:
	edit file test.nml

In the cuf version, there are two differents variants of the numerical scheme which are different in term of allocated memory. Use parameter implementationVersion = 0 or 1 to change.

Initial condition: a discontinuity along the domain diagonal

## Example of use

./euler2d_gpu ./test.nml

Output:
	VTK ascii file using VTK Image Data format

Visualization:
	paraview --data=euler2d_..vti

The cuf version was written in 2013 (with the PGI compiler), and the stdpar version in 2022.

In terms of performance, the stdpar version is able to run at more than 500 Mcell-update per seconds on A100, using large domain (2048x2048).
