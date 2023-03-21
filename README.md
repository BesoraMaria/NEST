NEST Volume NEST is a Python module designed to compute the volume occupied by the atoms of a molecule within a defined capsule shaped region. The occupancy of four quadrants Q1occ, Q2occ, Q3occ and Q4occ is also evaluated.
![General view](https://github.com/BesoraMaria/NEST/blob/main/general.png?raw=true)
The NEST tool is inspired on the great utility of SambVuca (Organometallics 2016, 35, 13, 2286–2293) to evaluate Buried Volumes and Molquo (Chem. Eur. J. 2012, 18, 995; Chem. Eur. J. 2012, 18, 14026) to evaluate quadrants occupancy. SambVca is designed to evaluate systems whose steric effects can be well explained within the volume of a sphere, we envisioned that a similar tool able to expand systems in one direction could be useful to evaluate the enantioselectivity of non-spherical (linear) systems. To do so, we define a volume with the shape of a capsule (cylinder with rounded ends), the length of the volume can be tuned by the user. The NEST is defined by a semi sphere (of 3.5 Å radii) plus a cylindrical unit of length defined by the user and ends with another semi sphere (of 3.5 Å radii). If the length is 0 the volume would be an sphere, for larger lengths the cylinder would be larger.
![Scheme](https://github.com/BesoraMaria/NEST/blob/main/figura.png?raw=true)
An atom “M” must be provided as origin of coordinates, usually the metal centre. Another atom, usually the central atom of the ligand “CL” must be provided to describe de y axis. And a third atom is necessary to define the z axis. This last atom will be only used in quadrant analysis. 
NEST needs of an xyz file, containing number of atoms in first line, and a line followed by cartessian coordinates (atom_type x y z). The results will depend on the atoms chosen to define the axis but also on the size of boxes (grid) and the scaling factor of Bondi Radii of atoms.

![Quadrants](https://github.com/BesoraMaria/NEST/blob/main/quatrecyan.png?raw=true)

NEST is available, as a code to be downloaded, but also is available as a web app in Streamlit: https://besoramaria-nest-nest-01-f2n5of.streamlit.app/

This progect has been developed by Dr. M Besora at the URV with the help of Dr. JI Mujika and G.D. Nunez.
More information is planned to be added soon
