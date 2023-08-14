# TwoTorsionSubgroups

This is the MAGMA code use the compute the 2-tosion subgroups in the arxiv preprint "Two-Torsion Subgroups of Some Modular Jacobians".

The paper presents a method to compute the 2-torsion subgroup of Jacobians of genus 3, 4 and 5 non-hyperelliptic curves. 


The code to compute the 2-torsion subgroup of the Jacobian of the genus 3 curve X0(75)/w25 is written in two files:
X075w25Bitangents.m - this computes the bitangents to the curve X0(75)/w25 
J075w25TwoTorsion.m - this uses bitangents computed above to compute the 2-torsion subgroup of the jacobian and its rational subgroup 
The file X075w25Bitangents.txt is also included, this simply states the bitangents computed using the first file.

The code used to compute J0(54)[2], the genus 4 example is written in two files:
J054TwoTorsionFOD.m - this computes a candidate for the field of definition of the 2-torsion subgroup, this is in fact the correct field of definition, as confirmed in the following file
J054TwoTorsion.m - this uses the field field of definition computed above to compute the 2-torsion subgroup and it's rational subgroup 
The file X054Tritangents.txt is also included, this simply stated the complete list of 120 tritangent planes to the curve. These are not explicitly used in these computations, but were used in a previous version of the paper, and it's simply included here for completness.

The files used for the genus 5 example J0(N), N =42, 55, 63, 72, 75 are as follows:
general.m - this file is used by all other X0NPlane.m files. It contains the general programs for the Newton-Raphson and lattice reduction steps. 
The file X0NQuadritangents.txt - simply states the orbits of quadritangents used in our computations 

X0NPlanes.m - computes some orbits of quadritangents 
J0NPlanes.m - uses the planes computed above to calculate the 2-torsion subgroup and it's rational subgroup

There is also a file X055Test.m - this demonstrates the fact that our computations of the minimal polynomials is more efficient than the MAGMA command MinimalPolynomial().
