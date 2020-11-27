# Kim & Hummer Model A (KH A)[1] & Hydrophobicity Scale Model (HPS)[2] LAMMPS Implementation
This project is an example of implementing continuous piecewise Lennard-Jones-like potentials in Lammps, such as KH A and HPS. 
Instead of writing a custom pair_style, only standard pair_styles are used. 
The primary downside of this method is that redundant pair_styles are required which lead to unnecessary inefficiency. 
However, this method can be used quickly and easily on most versions of LAMMPS without the need for being recompiled.
That stated, the authors have only tested this method on a limited number of LAMMPS versions and can't guarantee it will work with every version.

name|abc|tla|charge|radius|hps|type|mass
---|---|---|---|---|---|---|---
Alanine|A|ALA|0|5.04|0.73|1|71.08
Cysteine|C|CYS|0|5.48|0.595|2|103.1
Aspartic Acid|D|ASP|-1|5.58|0.378|3|115.1
Glutamic Acid|E|GLU|-1|5.92|0.459|4|129.1
Phenylalanine|F|PHE|0|6.36|1|5|147.2
Glycine|G|GLY|0|4.5|0.649|6|57.05
Histidine|H|HIS|0.5|6.08|0.514|7|137.1
Isoleucine|I|ILE|0|6.18|0.973|8|113.2
Lysine|K|LYS|1|6.36|0.514|9|128.2
Leucine|L|LEU|0|6.18|0.973|10|113.2
Methionine|M|MET|0|6.18|0.838|11|131.2
Asparagine|N|ASN|0|5.68|0.432|12|114.1
Proline|P|PRO|0|5.56|1|13|97.12
Glutamine|Q|GLN|0|6.02|0.514|14|128.1
Arginine|R|ARG|1|6.56|0|15|156.2
Serine|S|SER|0|5.18|0.595|16|87.08
Threonine|T|THR|0|5.62|0.676|17|101.1
Valine|V|VAL|0|5.86|0.892|18|99.07
Tryptophan|W|TRP|0|6.78|0.946|19|186.2
Tyrosine|Y|TYR|0|6.46|0.865|20|163.2

1. Kim YC, Hummer G. (2008). Coarse-grained models for simulation of multiprotein complexes: application to ubiquitin binding. J. Mol. Biol. 375,5 1416-1433. https://doi.org/10.1016/j.jmb.2007.11.063
2. G.L. Dignon, W. Zheng, Y.C. Kim, R.B. Best, J. Mittal. (2018). Sequence determinants of protein phase behavior from a coarse-grained model. PLoS Comput. Biol. https://doi.org/10.1371/journal.pcbi.1005941
