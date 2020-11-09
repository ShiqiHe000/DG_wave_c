Adaptive Mesh Refinement and Coarsening
************************************************
Modern supercomputers enable scientists and engineers to 
tackle large scale problems. 
However, applying a static uniform fine mesh to the whole computational domain 
can be extremely expensive for a grand challenge problem, and sometimes impractical.

Further, many CFD problems aim to track features that 
are much smaller than the overall scale of the problem. 
In some problems high resolutions are only required for computationally difficult regions.

An intuitive strategy to tackle this problem is to generate 
high-resolution meshes in needed computational regions. 
This strategy is the so-called **adaptive mesh refinement and coarsening (AMR)**. 


Adptivity: hp-refinement
-------------------------------------------------
.. list-table:: Two types of refinement

    * - .. figure:: /image/AMR/h_refinement.png 
	
	   h-refinement		

      - .. figure:: /image/AMR/p_refinement.png

	   p-refinement

In this work, two types of refinement methods are invoked: 
**h-refinement** and **p-refinement**. 
h-refinement subdivides a target element into smaller children elements.
p-type refinement is exclusive to high-order methods. 
The mesh resolution can be improved by increasing its polynomial order with respect to the target direction.

Decisions between h- and p-type refinement
```````````````````````````````````````````````````
A strategy to take advantage from both of the refinement methods: 
when the solution is poor, 
h-refinement should be imposed to eliminate the numerical error. 
Once a smooth solution is obtained, 
exponential convergence can be performed by prescribing p-refinement.

- Poor solution

  * h-refinement

- Smooth solution

  * p-refinement


Error Estimator
-----------------------------------
How to evaluate the smoothness of the solution? 

- Use an error estimator!

An *a posteriori* error estimator :cite:`Mavriplis_refinement` is implemented in this work. 



