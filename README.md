# Implementation of Kim & Hummer (KH) Model [1] in LAMMPS 
This project uses continuous piecewise Lennard-Jones-like potentials to implement the KH Model A in LAMMPS. 
Instead of writing a custom pair_style, only standard pair_styles are used. 
This method can thus be used quickly and easily on most versions of LAMMPS, which however may lead to reduction of efficiency compared to compiled custom potentials. 
That stated, the authors have only tested this method on a limited number of LAMMPS versions and can't guarantee it will work with every version.
Provided python scripts can be used to initialize simulations from pdb files. 

The code of this project was used in the simulations in [2], which can be used a reference. 

Table of important parameters where listed type is indentifier used in LAMMPS input files:

name|abc|tla|charge (e)|radius (Å)|type|mass (g/mole)
---|---|---|---|---|---|---
Alanine|A|ALA|0|5.04|1|71.08
Cysteine|C|CYS|0|5.48|2|103.1
Aspartic Acid|D|ASP|-1|5.58|3|115.1
Glutamic Acid|E|GLU|-1|5.92|4|129.1
Phenylalanine|F|PHE|0|6.36|5|147.2
Glycine|G|GLY|0|4.5|6|57.05
Histidine|H|HIS|0.5|6.08|7|137.1
Isoleucine|I|ILE|0|6.18|8|113.2
Lysine|K|LYS|1|6.36|9|128.2
Leucine|L|LEU|0|6.18|10|113.2
Methionine|M|MET|0|6.18|11|131.2
Asparagine|N|ASN|0|5.68|12|114.1
Proline|P|PRO|0|5.56|13|97.12
Glutamine|Q|GLN|0|6.02|14|128.1
Arginine|R|ARG|1|6.56|15|156.2
Serine|S|SER|0|5.18|16|87.08
Threonine|T|THR|0|5.62|17|101.1
Valine|V|VAL|0|5.86|18|99.07
Tryptophan|W|TRP|0|6.78|19|186.2
Tyrosine|Y|TYR|0|6.46|20|163.2

Miyazawa-Jernigan contact potential between residues i & j in RT [3]:

0|1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20
---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
 1|-2.72|-3.57|-1.70|-1.51|-4.81|-2.31|-2.41|-4.58|-1.31|-4.91|-3.94|-1.84|-2.03|-1.89|-1.83|-2.01|-2.32|-4.04|-3.82|-3.36
 2|     |-5.44|-2.41|-2.27|-5.80|-3.16|-3.60|-5.50|-1.95|-5.83|-4.99|-2.59|-3.07|-2.85|-2.57|-2.86|-3.11|-4.96|-4.95|-4.16
 3|     |     |-1.21|-1.02|-3.48|-1.59|-2.32|-3.17|-1.68|-3.40|-2.57|-1.68|-1.33|-1.46|-2.29|-1.63|-1.80|-2.48|-2.84|-2.76
 4|     |     |     |-0.91|-3.56|-1.22|-2.15|-3.27|-1.80|-3.59|-2.89|-1.51|-1.26|-1.42|-2.27|-1.48|-1.74|-2.67|-2.99|-2.79
 5|     |     |     |     |-7.26|-4.13|-4.77|-6.84|-3.36|-7.28|-6.56|-3.75|-4.25|-4.10|-3.98|-4.02|-4.28|-6.29|-6.16|-5.66
 6|     |     |     |     |     |-2.24|-2.15|-3.78|-1.15|-4.16|-3.39|-1.74|-1.87|-1.66|-1.72|-1.82|-2.08|-3.38|-3.42|-3.01
 7|     |     |     |     |     |     |-3.05|-4.14|-1.35|-4.54|-3.98|-2.08|-2.25|-1.98|-2.16|-2.11|-2.42|-3.58|-3.98|-3.52
 8|     |     |     |     |     |     |     |-6.54|-3.01|-7.04|-6.02|-3.24|-3.76|-3.67|-3.63|-3.52|-4.03|-6.05|-5.78|-5.25
 9|     |     |     |     |     |     |     |     |-0.12|-3.37|-2.48|-1.21|-0.97|-1.29|-0.59|-1.05|-1.31|-2.49|-2.69|-2.60
10|     |     |     |     |     |     |     |     |     |-7.37|-6.41|-3.74|-4.20|-4.04|-4.03|-3.92|-4.34|-6.48|-6.14|-5.67
11|     |     |     |     |     |     |     |     |     |     |-5.46|-2.95|-3.45|-3.30|-3.12|-3.03|-3.51|-5.32|-5.55|-4.91 
12|     |     |     |     |     |     |     |     |     |     |     |-1.68|-1.53|-1.71|-1.64|-1.58|-1.88|-2.83|-3.07|-2.76 
13|     |     |     |     |     |     |     |     |     |     |     |     |-1.75|-1.73|-1.70|-1.57|-1.90|-3.32|-3.73|-3.19
14|     |     |     |     |     |     |     |     |     |     |     |     |     |-1.54|-1.80|-1.49|-1.90|-3.07|-3.11|-2.97
15|     |     |     |     |     |     |     |     |     |     |     |     |     |     |-1.55|-1.62|-1.90|-3.07|-3.41|-3.16
16|     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |-1.67|-1.96|-3.05|-2.99|-2.78
17|     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |-2.12|-3.46|-3.22|-3.01
18|     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |-5.52|-5.18|-4.62
19|     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |-5.06|-4.66 
20|     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |-4.17

1. Kim YC, Hummer G. (2008). Coarse-grained models for simulation of multiprotein complexes: application to ubiquitin binding. J. Mol. Biol. 375,5 1416-1433. https://doi.org/10.1016/j.jmb.2007.11.063
2. H. Smith, N. Pinkerton, D. B. Heisler, E. Kudryashova, A. R. Hall, K. R. Karch, A. Norris, V. Wysocki, M. Sotomayor, E. Reisler, D. Vavylonis, D. S. Kudryashov (2020) Rounding out the understanding of ACD toxicity with the discovery of cyclic forms of actin oligomers. Int. J. Mol. Sci. 2021, 22, 718; https://doi.org/10.3390/ijms22020718
3. Miyazawa, S., Jernigan, R.L. (1996). Residue–Residue Potentials with a Favorable Contact Pair Term and an Unfavorable High Packing Density Term, for Simulation and Threading. J. Mol. Biol., 256,3 623-644. https://doi.org/10.1006/jmbi.1996.0114

