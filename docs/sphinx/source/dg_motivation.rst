Motivation
============================
.. figure:: /image/motivation/airplane_turbulence.jpg
	:alt: "Turbulent flow around an airplane." 
	:align: center	

	Turbulent flow around an airplane. 

Scientists and engineers want to study **complex fluid flow**, *e.g.*, turbulence, 
which can help them to design more fuel-efficient airplane. 

To study complex fluid flows, you need a **high-order** (high accuracy) method to capture the details inside fluid flow. 
In our work, we use the 

- **discontinuous spetral element method (DG-SEM)**. 

One potent method to fully resolve turbulent flow is to use **direct numerical simulation (DNS)**. 
However, DNS is very computational expensive: the computational resources required by a DNS would exceed the capacity of the most powerful computers currently available. 
To save the computational cost we only bring high resolution to the computational difficult regions to the mesh:

- **parallel adaptive mesh refinement (AMR)** is appltied to our solver. 

.. figure:: /image/motivation/amr_airplane.png
	:align: center
	:figwidth: 70%

	AMR example:iosbars and mesh cut on a business jet configuration computed with AMR approach.

However, parallel AMR imposed a huge challange on data encoding and managing. 
We propose to use **hash table** data structure for AMR data management 
to supplant traditional tree-structured AMR. 

.. figure:: /image/motivation/hash_vs_tree.png
	:align: center
	:figwidth: 70%

	Tree-structured AMR and hash table AMR features comparison.

Another challange is caused by the adptivity of the solver.
Dynamic mesh introduces load imbalance among the processors: 
processors with heavier workload need more time to performs computation than the ones with less workload. 
Load imbalance can largely degrade the performance of supercomputers. 
To tackle this problem, a **space-filling curve (SFC) based repartitioning algorithm**
is applied to redistribute the workload among processors evenly and periodicially. 
The SFC we choose is called **Hilbert curve**. 

.. figure:: /image/motivation/Hilbert_different_level.png
	:align: center
	:figwidth: 30%

	The Hilbert curve traverses the mesh with mixing levels.
