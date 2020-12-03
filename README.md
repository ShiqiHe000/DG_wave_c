# Dynamic Load Balancing for a hp-adaptive Discontinuous Galerkin Wave Equation Solver via Space-Filling Curve and Advanced Data Structure

[![Documentation Status](https://readthedocs.org/projects/dg-wave-c/badge/?version=latest)](https://dg-wave-c.readthedocs.io/en/latest/?badge=latest)

<!--ts-->
   * [Dynamic Load Balancing for a hp-adaptive Discontinuous Galerkin Wave Equation Solver via Space-Filling Curve and Advanced Data Structure](#dynamic-load-balancing-for-a-hp-adaptive-discontinuous-galerkin-wave-equation-solver-via-space-filling-curve-and-advanced-data-structure)
   * [Table of Contents](#table-of-contents)
      * [Introduction](#introduction)
      * [Setup](#setup)
      * [Documentation](#documentation)
      * [Source Code Documentation](#source-code-documentation)
      * [Approximation of Wave Equation](#approximation-of-wave-equation)

<!-- Added by: shiqi, at: Wed Dec  2 16:55:06 EST 2020 -->

<!--te-->

## Introduction
We combine a high-order method -- the **discontinuous Galerkin [spectral element method](https://en.wikipedia.org/wiki/Spectral_element_method) (DG-SEM)**, 
with parallel [**adaptive mesh refinement and coarsening (AMR)**](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement) techniques and apply it to a **two-dimensional [wave equation](https://en.wikipedia.org/wiki/Wave_equation) solver**.

Advanced data structures and dynamic load balancing are applied to the solver to achieve efficient data management and high-level parallelism. 

## Setup
### What You Need
* C++
* [CMake](https://cmake.org/) (at least version 3.9)
* [GCC](https://gcc.gnu.org/) 7.5.0 (GNU Compiler Collection)
* [OpenMPI](https://www.open-mpi.org/) 4.0.2

### Compile and Execute
```
mkdir build && cd build
cmake ..
make 
mpirun -np 4 main
```
You also can use parallel `make`, *e.g.*, `make -j 16`. 16 threads will be used to compile the code. 

`mpirun -np 4 main` executes the program with 4 processors, you could change the processor number as you want.

## Documentation
A detailed documentation an be found at [here](https://dg-wave-c.readthedocs.io/en/latest/).

## Source Code Documentation
Source code explanation:
[Source code documentation]( https://shiqihe000.github.io/DG_wave_c/doxygen/html/index.html)

## Approximation of Wave Equation
The basic model of wave propagation is the wave equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" title="\frac{\partial ^{2}p}{\partial t^{2}}-c^{2}(p_{xx}+p_{yy})=0" /></a>

The variable `p` represents the acoustic pressure and `c` is the sound speed. 

## AMR: hp-adaptivity
Two types of refinements are implemented in this work: h-refinement and p-refinement. 
### h-refinement
<p align="center">
  <img src="./imgs/h_refinement.png" width="100" height = "40" >
</p>
Subdivide an element into children elements. 
### p-refinement
<p align="center">
  <img src="./imgs/p_refinement.png" width="100" height = "40" >
</p>
Raise polynomial orders inside the targeted element. 
