# Kim & Hummer Model A (KH A)[1] & Hydrophobicity Scale Model (HPS)[2] LAMMPS Implementation
This project is an example of implementing continuous piecewise Lennard-Jones-like potentials in Lammps, such as KH A and HPS. 
Instead of writing a custom pair_style, only standard pair_styles are used. 
The primary downside of this method is that redundant pair_styles are required which lead to unnecessary inefficiency. 
However, this method can be used quickly and easily on most versions of LAMMPS without the need for being recompiled.
That stated, the authors have only tested this method on a limited number of LAMMPS versions and can't guarantee it will work with every version.

1. Kim YC, Hummer G. (2008). Coarse-grained models for simulation of multiprotein complexes: application to ubiquitin binding. J. Mol. Biol. 375,5 1416-1433. https://doi.org/10.1016/j.jmb.2007.11.063
2. G.L. Dignon, W. Zheng, Y.C. Kim, R.B. Best, J. Mittal. (2018). Sequence determinants of protein phase behavior from a coarse-grained model. PLoS Comput. Biol. https://doi.org/10.1371/journal.pcbi.1005941
