# Discontinuous Galerkin Wave Equation Slover

[![Build Status](https://travis-ci.com/ShiqiHe000/DG_wave_c.svg?branch=master)](https://travis-ci.com/ShiqiHe000/DG_wave_c)

[![Documentation Status](https://readthedocs.org/projects/dg-wave-c/badge/?version=latest)](https://dg-wave-c.readthedocs.io/en/latest/?badge=latest)

## Documentation
A detailed documentation an be found at [https://2d-advection.readthedocs.io/en/latest/]

## Documentation from the source code
[Source code documentation]( https://shiqihe000.github.io/DG_wave_c/doxygen/html/index.html)

## Approximation of wave equation
The basic model of wave propgation is the wave equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;^{2}p}{\partial&space;t^{2}}-c^{2}(p_{xx}&plus;p_{yy})=0" title="\frac{\partial ^{2}p}{\partial t^{2}}-c^{2}(p_{xx}+p_{yy})=0" /></a>

The variable `p` might represent the acoustic pressure and `c` would be the sound speed. 
