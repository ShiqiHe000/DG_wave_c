Spectral Approximation on a square
*********************************************

Approximation of Wave Propagation
=============================================

Basic Model
---------------------------------------------
The basic model is the linear wave equation with the form:

.. math::
        \frac{\partial ^2 p}{\partial t^2} - c^2 (\frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2}) = 0

The wave equation is the fundamental equation of acoustics. 
It is based on two improtant approximation, namely, that the flow may be treated as *inviscide* and that *convective derivatives are negligible in comparison to unsteady derivatives*.
(we neglect viscous and other diffusion effect(heat), when convection transfer is much faster than diffusion transfer of mass, momentum or energy.)

The variable :math:`p` may represent acoustic pressure in an otherwise quiescent gas and :math:`c` could be sound speed. 

In order to solve the second order equation, we re-write the equation as a system of three first order equations.

Convert the wave equation to a system of first order equation, let:

.. math::
        u_t = - p_x,v_t = -p_y.

:math:`u` and :math:`v` correspond to the components of the velocity in a fluid flow. 

Assuming the order of mixed partial derivatives does not matter, then:

.. math::
        \frac{\partial^2 p}{\partial t^2} + c^2((u_x)_t + (v_y)_t) = 0.

Combining with initial conditions,

.. math::
        p_t + c^2(u_x + v_y) = 0.

We now obtain the system of equations by grouping the equation for pressure and two velocity components

.. math::

        \begin{bmatrix}
        p\\ 
        u\\ 
        v
        \end{bmatrix}_t +
        \begin{bmatrix}
        0& c^2 & 0\\ 
        1& 0 & 0\\ 
        0& 0 & 0
        \end{bmatrix}
        \begin{bmatrix}
        p\\ 
        u\\ 
        v
        \end{bmatrix}_x+
        \begin{bmatrix}
        0 & 0 & c^2\\ 
        0& 0 & 0\\ 
        1&  0& 0
        \end{bmatrix}\begin{bmatrix}
        p\\ 
        u\\ 
        v
        \end{bmatrix}_y
 
or 

.. math::
        \mathbf{q_t} + A\mathbf{q_x} +B\mathbf{q_y} = 0

Since :math:`A` and :math:`B` are constants, we can bring them inside the derivatives

.. math::
        \mathbf{q_t} + \mathbf{f_x} + \mathbf{g_y} = 0 \\
        \mathbf{f_x} = A\mathbf{q_x} \\
        \mathbf{g_y} = B\mathbf{q_y} \\

This is known as **Conservation law** form since it can be written as 

.. math::
        \mathbf{q_t} + \bigtriangledown \cdot F = 0

where the vector flux :math:`F = \mathbf{f}\widehat{x}+\mathbf{g}\widehat{y}`. 

The term conservation law follows from the fact that the differential equation is what we get when we apply the divergence theorem to the integral conservation law.

.. math::
        \frac{d}{dt} \int_{V} \mathbf{q}dV = - \int_{S} F \cdot \widehat{n} dS

Riemann Problem for Conservation Law
---------------------------------------------

Introduction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A `Riemann problem`_, named after Bernhard Riemann, is a specific initial value problem composed of a conservation equation together with piecewise constant initial data which has a single discontinuity in the domain of interest. The Riemann problem is very useful for the understanding of equations like Euler conservation equations because all properties, such as shocks and rarefaction waves, appear as characteristics in the solution. It also gives an exact solution to some complex nonlinear equations, such as the Euler equations. 

.. _`Riemann problem`: https://en.wikipedia.org/wiki/Riemann_problem

.. image:: /image/Riemann1.png





Riemann Solver
^^^^^^^^^^^^^^^^^^^^^^^^
Here we build a Riemann problem for the hyperbolic, constant coefficient system with proper initial condition. 

.. math::
        \mathbf{q_t} + A\mathbf{q_x} +B\mathbf{q_y} = 0

The coefficient matrices :math:`A` and :math:`B` have :math:`m` real eigenvalues :math:`\lambda_i` and :math:`m` linearly independent eigenvectors :math:`\mathbf{K}^{(i)}`, where :math:`m` is the equation number :doc:`./dg_reference`.  






The Nodal Discontinuous Galerkin Approximation
-----------------------------------------------
We will implement the discontinuous Galerkin spectral element approximation of two-dimensional conservation law on a square domain.

.. math::
        \mathbf{q_t} + \mathbf{f_x} +\mathbf{g_y}= 0, x \in (L, R), y \in(D, U) 
        :label: equ1
        
The spectral element approximation starts with a weak form of :eq:`equ1`. We multiply :eq:`equ1` by a test function, integrate and subdivide into elements

.. math::
        \sum_{k=1}^{K}\left [ \int_{x_{k-1}}^{x_k} (\mathbf{q}_t+\mathbf{f}_x + \mathbf{g}_y)\phi dx\right ] = 0
        :label: equ2

We map :eq:`equ2` onto reference space by **affine map** :eq:`equ3`

.. math::
        x = x_{k-1} + \frac{\xi +1}{2} \Delta x_k, \Delta x_k = x_k - x_{k+1}\\
        y = y_{k-1} + \frac{\eta  +1}{2} \Delta y_k, \Delta y_k = y_k - y_{k+1}\\
        dx = \frac{\Delta x_k}{2}d\xi , \frac{\partial}{\partial x} = \frac{2}{\Delta x_k}\frac{\partial }{\xi}
        :label: equ3

The solution and fluxes are approximated by polynomials of degree N and represent the polynomials in nodal, Lagrange form

.. math::
        \mathbf{q} \approx \mathbf{Q} = \sum_{n=0}^{N}\sum_{m=0}^{M}\mathbf{Q}_{n,m} l_n(x)l_m(y)\\
        \mathbf{F}_{n,m}\widehat{x} + \mathbf{G}_{n,m}\widehat{y} = B\mathbf{Q}_{n,m}\widehat{x} + C\mathbf{Q}_{n,m}\widehat{y}
        :label: equ4

where :math:`\mathbf{F}_{n,m}\widehat{x} + \mathbf{G}_{n,m}\widehat{y} = B\mathbf{Q}_{n,m}\widehat{x} + C\mathbf{Q}_{n,m}\widehat{y}`. We subsitute the approximations into the weak form of the PDE, and let :math:`(\mathbf{Q}_t, \phi _{ij}) + (\bigtriangledown \cdot \mathbf{F}, \phi_{ij}) = 0.`

If we apply  Green's identity to the second intergal

.. math::
        (\bigtriangledown \cdot \mathbf{F}, \phi_{ij}) = \int_{l}^{r} \phi_{ij} \bigtriangledown \cdot \mathbf{F}dxdy = \frac{\Delta x}{2} \int_{-1}^{1}\phi_{ij} \mathbf{f}_{\xi } d \xi +  \frac{\Delta y}{2} \int_{-1}^{1}\phi_{ij} \mathbf{g}_{\eta } d \eta 

The Nurmerical flux
----------------------------------------------

Time Integration
-----------------------------------------------

Change of Interval
---------------------------------------------

Benchmark Solution: Plane wave Propagation
----------------------------------------------
We represent a plane Gaussian wave through the grid. 

The plane wave is defined as:

.. math::
        \begin{bmatrix}
        p\\ 
        u\\ 
        v
        \end{bmatrix} =
        \begin{bmatrix}
        1\\ 
        \frac{k_x}{c}\\ 
        \frac{k_y}{c}
        \end{bmatrix}
        e^{-\frac{(k_x(x-x_0)+k_y(y-y_0)-ct)^2}{d^2}}

Where :math:`\mathbf{k}` is the wavevector and it is normalized to satisfiey :math:`k_x^2 + k_y^2 = 1`.
The wavevector is choosen as :math:`\mathbf{k} = (\sqrt{2}/2, \sqrt{2}/2)`
This is a wave with Gaussian shape where we compute the parameter :math:`d` from the full width at half maximum, :math:`\omega  = 0.2`, by math:`d = \omega/2\sqrt{ln2}`. 
The other parameters are :math:`c = 1` and :math:`x_0 = y_0 = -0.8`. 

Performance Evaluation
-------------------------------------------
Exact boundary solutions are imposed on the 4 side of the computation domain. The initial condition is setting `t=0.0` of the exact solution. 

1 element 
^^^^^^^^^^^^^^^^^^^^^^
Domain: :math:`x \in [0.0, 1.0], y\in [0.0, 1.0]`.

Time step: :math:`\Delta t = 2.0\times 10^{-4}`

Fig(1), shows the error performances. 

.. image:: /image/2d_1_element_error.png

4 element2
^^^^^^^^^^^^^^^^^^^^^^^
Domain: :math:`x \in [0.0, 1.0], y\in [0.0, 1.0]`.

Time step: :math:`\Delta t = 2.0\times 10^{-4}`

Fig(2), shows the error performances.

.. image:: /image/2d_4_elements.png

16 elements
^^^^^^^^^^^^^^^^^^^
Domain: :math:`x \in [0.0, 1.0], y\in [0.0, 1.0]`.

Time step: :math:`\Delta t = 1.0\times 10^{-5}`

Fig(3), shows the error performances.

.. image:: /image/2d_16_elements_error.png

64 elements
^^^^^^^^^^^^^^^^^^^^^^^
Domain: :math:`x \in [0.0, 8.0], y\in [0.0, 8.0]`.

Time step: :math:`\Delta t = 1.0\times 10^{-5}`

Fig(3), shows the error performances.

.. image:: /image/2d_64_elements_error.png
