![euler2d 350x350](https://github.com/pkestene/euler2d_cudaFortran/blob/master/euler2d.png?raw=true)

# euler2d_cudaFortran

2nd order Godunov solver for 2d Euler equations written in CUDA Fortran

Note: If your more interested into c++ for GPUs, have a look at version using C++/kokkos for performance portability:
https://github.com/pkestene/euler2d_kokkos

## Short description

This code solves the 2D Euler equations in a regular cartesian mesh
using a 2nd order godunov-based finite volume scheme.

There are actually two version, one is using the so called [cuda fortran](https://developer.nvidia.com/cuda-fortran) programming model, and the other is designed with the newer [fortran standard parallelism](https://developer.nvidia.com/blog/accelerating-fortran-do-concurrent-with-gpus-and-the-nvidia-hpc-sdk/) (stdpar).

## Nvhpc compiler for Cuda Fortran and stdpar

You need to have installed [nvhpc] to build this application. then, e.g

```shell
module load nvhpc/23.11
```

By default, nvfortran will detect visible GPU available on the system, and generate code for it.

```shell
# for the cuda fortran version
cd cuf
make
# for the stdpar version
cd stdpar
make
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

In terms of performance,
- stdpar version is able to run at more than 500 Mcell-updates per seconds on A100, using large domain (2048x2048).
- cuda-fortran version has roughly the same performance (maybe slightly less performant).
