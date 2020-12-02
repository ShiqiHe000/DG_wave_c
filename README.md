# Dynamic Load Balancing for a hp-adaptive Discontinuous Galerkin Wave Equation Solver via Space-Filling Curve and Advanced Data Structure

[![Documentation Status](https://readthedocs.org/projects/dg-wave-c/badge/?version=latest)](https://dg-wave-c.readthedocs.io/en/latest/?badge=latest)

<!--ts-->
Table of Contents
=================

   * [Dynamic Load Balancing for a hp-adaptive Discontinuous Galerkin Wave Equation Solver via Space-Filling Curve and Advanced Data Structure](#dynamic-load-balancing-for-a-hp-adaptive-discontinuous-galerkin-wave-equation-solver-via-space-filling-curve-and-advanced-data-structure)
      * [Introduction](#introduction)
      * [Documentation](#documentation)
      * [Source Code Documentation](#source-code-documentation)
      * [Approximation of Wave Equation](#approximation-of-wave-equation)
<!--te-->

## Introduction
We combine a high-order method -- the **discontinuous Galerkin [spectral element method](https://en.wikipedia.org/wiki/Spectral_element_method) (DG-SEM)**, 
with parallel [**adaptive mesh refinement and coarsening (AMR)**](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement) techniques and apply it to a **two-dimensional [wave equation](https://en.wikipedia.org/wiki/Wave_equation) solver**.

Advanced data structures and dynamic load balancing are applied to the solver to achieve efficient data management and high-level parallelism. 

## Documentation
A detailed documentation an be found at [here](https://dg-wave-c.readthedocs.io/en/latest/).

## Source Code Documentation
[Source code documentation]( https://shiqihe000.github.io/DG_wave_c/doxygen/html/index.html)

## Approximation of Wave Equation
The basic model of wave propagation is the wave equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" title="\frac{\partial ^{2}p}{\partial t^{2}}-c^{2}(p_{xx}+p_{yy})=0" /></a>

The variable `p` might represent the acoustic pressure and `c` would be the sound speed. 
