## Introduction
NEST is a Python module designed to compute the volume occupied by the atoms of a molecule within a defined capsule shaped region. The occupancy of four quadrants Q1occ, Q2occ, Q3occ and Q4occ is also evaluated.
![General view](https://github.com/BesoraMaria/NEST/blob/main/general.png?raw=true)
The NEST tool is inspired on the great utility of SambVuca (Organometallics 2016, 35, 13, 2286–2293) to evaluate Buried Volumes and Molquo (Chem. Eur. J. 2012, 18, 995; Chem. Eur. J. 2012, 18, 14026) to evaluate quadrants occupancy. SambVca is designed to evaluate systems whose steric effects can be well explained within the volume of a sphere, we envisioned that a similar tool able to expand systems in one direction could be useful to evaluate the enantioselectivity of non-spherical (linear) systems. To do so, we define a volume with the shape of a capsule (cylinder with rounded ends), the length of the volume can be tuned by the user. The NEST is defined by a semi sphere (of 3.5 Å radii) plus a cylindrical unit of length defined by the user and ends with another semi sphere (of 3.5 Å radii). If the length is 0 the volume would be an sphere, for larger lengths the cylinder would be larger.
![Scheme](https://github.com/BesoraMaria/NEST/blob/main/figure.png?raw=true)
An atom “M” must be provided as origin of coordinates, usually the metal centre. Another atom, usually the central atom of the ligand “CL” must be provided to describe de y axis. And a third atom is necessary to define the z axis. This last atom will be only used in quadrant analysis. 
NEST needs of an xyz file, containing number of atoms in first line, and a line followed by cartessian coordinates (Atomic_symbol x y z). The results will depend on the atoms chosen to define the axis but also on the size of boxes (grid) and the scaling factor of Bondi Radii of atoms.
                   ![Quadrants](https://github.com/BesoraMaria/NEST/blob/main/quatrecyan.png?raw=true)

The volume occupation for each of the 4 quadrants, defined from the metal center until the end of the capsule (oposite direction to CL) and taking the atom defined as Z to set position of quadrant boundaries is also computed. Quadrant occupation is given as Q1occ, Q2occ, Q3occ and Q4occ.

## Availability
NEST is available, as a code to be downloaded, but also is available as a web app in Streamlit: https://besoramaria-nest-nest-01-f2n5of.streamlit.app/

## Citing NEST
G. Zuccarello, L.J. Nannini, A. Arroyo-Bondía, N. Fincias, I. Arranz, A.H. Pérez-Jimeno, M. Peeters, I. Martín-Torres, A. Sadurní, V. García, Y. Wang, M.S. Kirillova, M. Montesinos-Magraner, U. Caniparoli, G.D. Núñez, F. Maseras, M. Besora, I. Escofet, A.M. Echavarren. Enantioselective Catalysis with Pyrrolidinyl Gold(I) Complexes: DFT and NEST Analysis of the Chiral Binding Pocket. * *JACS Au* *, 2023, XXXX. https://doi.org/10.1021/jacsau.3c00159

## Usage
To use NEST App the user needs an xyz file of the molecules of interest. The xyz file should have the standard format of first line giving the total number of atoms, second line with title or empty line, and coordinates starting from the third line. Coordinates should be expressed as:

atom_symbol   x_coordinate   y_coordinate   z_coordinate

Then fill in the atom number to place the origin (M), atom number of the ligand to place the y axis (CL) and atom number (Z) to place the z axis. This is a crucial step, if the orientation is not the desired, the description obtained will not have the expected meaning. In more detail:
 
The NEST App will ask you "Number of origin atom (Metal)"
 
The atom selected as "origin atom" will become the origin of coordinates (0,0  0,0 0,0) after running the app. This atom is expected to be, usually, the considered metal center, although it can be any atom. The number corresponds to the number of the selected atom in the atom list (first atom 1, second atom 2, and so on). A ghost atom could be included to place the center of coordinates at the desired point. The atom number can be easily checked by any visualizer able to open xyz files, like jmol, vmd, avogadro,...

The NEST App will ask you "Number of atom in y axis (Central Ligand):"
 
This "atom in y axis” together with the origin atom given above, will give the y axis. This atom will be on the y negative axis, and hence the coordinates of this atom after running the program will be: (0,0  -distance 0,0) where the distance is the distance between this atom and the atom marked as origin. This atom is expected to be, usually, the central atom of the ligand coordinated to the metal center, but this is not mandatory, and it can be any atom. The number corresponds to the number of the selected atom in the atom list (first atom 1, second atom 2, and so on). A ghost atom could be included to place the center of coordinates at the desired point. The atom number can be easily checked by any visualizer able to open xyz files, like jmol, vmd, avogadro,...

The NEST App will ask you "Number of atom in z axis (out of xy plane):" This "atom in z axis” together with the xy plane defined above, will give the z axis. This atom will be on the z negative axis, and hence the coordinates of this atom after running the program will be: (0,0  value1  value2). This atom is expected to be, usually, an atom of the ligand close to be perpendicular to the xy plane, but this is not mandatory, and it can be any atom. The number corresponds to the number of the selected atom in the atom list (first atom 1, second atom 2, and so on). A ghost atom could be included to place the center of coordinates at the desired point. The atom number can be easily checked by any visualizer able to open xyz files, like jmol, vmd, avogadro,...

Also the NEST App will ask you "Atoms excluded in V calculation (enter atom number):" If you need any atom or atoms to be excluded in the volume calculation, their atom numbers must be introduced in the format, number space number. A single atom would work (i.e. 86 ) as well as many (i.e. 86 87 88 89 ). Usually we would delete atoms present in the structure that are not relevant in our study, such as some small ligand, solvent molecules or other. In case of added ghost atoms they should be deleted in this section.

Now is time to set the lenght of the NEST (default=2): The value you introduce here is the value of the "lenght" of the cylinder, the overall lenght of the nest is this value plus 7 Å. For a value of zero, you get an sphere of 3.5 radius (hence results are equal within accuracy to those of SambVca) for values of 3 Å, we would get a NEST of 10 Å lenght. Measuring the distances of the regions of interest to the origin is important to set this length correctly. 

The atomic radii are important for these calculations, depending on them there will be more or less occupied volume. Based in SambVca descriptor, the NEST app uses Bondi radii, by default they are scales by 1.17 (as used in SambVca, see bibliography),but scaling factor can be tuned depending on system needs.

Finally, the user needs to define the grid of the calculation of the boxes, note that the desire value is 0.1, and can go to 0.5 (calculation faster, but losing precision), if the value is smaller than 0.1, the calculation will bcan become very demanding.

## Results
Once the calculations finished the user will obtaing 5 files and 5 numerical results:

The five numerical values are the NEST occupied volume (all capsule) and the volume occupied in four different 3D quadrants (only from M to the end of the capsule) Q1occ, Q2occ, Q3occ and Q4occ. 
Also, five different .xyz files can be downloaded, one with the full NEST volume represented. Please note that points represent space not occupied by atoms. And four with one 3D representation for each quadrant. For this files we recommend to use the Jmol interface to visualize them. Please note that they can be very large files if grid is small and length large. In more detail:

NEST occupied volume in %, it gives the occupied area respect to the total area (occupied and free). This data is printed on screen. The free volume of the NEST can be easily observed as each empty point will be represented in the generated new file called geomNEST_namefile.xyz . This file can be opened with jmol, VMD or similar. Please note that if the used grid is small this file will contain a large number of points and hence some programs will have serious problems to open it (Avogadro, Molden, ..). We strongly recommend opening this file with jmol to check that the orientation was correctly taken and the atoms to be removed were correctly selected. In the title section of the xyz file, the NEST occupied volume, and the Q1 to Q4 occupied volumes in % are reported.

Q1occ(x+,y+,z+) in %, it gives the occupied area respect to the  Q1 area (occupied and free). This data will be printed on the screen. Area Q1 is defined by positive x, y and z. Is just the quadrant infront of the CL-M and it does not include the parts of the ligand or molecule that are behind M (origin of atoms). The Q1 area and the free/occupied volume can be checked by inspection of the geomNEST_Q1_namefile.xyz file generated  with jmol, VMD or similar.
  
Q2occ(x-,y+,z+)  in %, it gives the occupied area respect to the  Q2 area (occupied and free). This data will be printed on the screen. Area Q2 is defined by positive y and z and negative x. Is just the quadrant infront of the CL-M and it does not include the parts of the ligand or molecule that are behind M (origin of atoms). The Q2 area and the free/occupied volume can be checked by inspection of the geomNEST_Q2_namefile.xyz file generated with jmol, VMD or similar.

Q3occ(x+,y+,z-) in %, it gives the occupied area respect to the  Q3 area (occupied and free). This data will be printed on the screen. Area Q3 is defined by positive x and y and negative z. Is just the quadrant in front of the CL-M and it does not include the parts of the ligand or molecule that are behind M (origin of atoms). The Q3 area and the free/occupied volume can be checked by inspection of the file generated geomNEST_Q3_namefile.xyz with jmol, VMD or similar.

Q4occ(x-,y+,z-) in %, it gives the occupied area respect to the  Q4 area (occupied and free). This data will be printed on the screen. Area Q4 is defined by positive y and negative x and z. Is just the quadrant infront of the CL-M and it does not include the parts of the ligand or molecule that are behind M (origin of atoms). The Q4 area and the free/occupied volume can be checked by inspection of the geomNEST_Q4_namefile.xyz file generated with jmol, VMD or similar

## Development

The NEST App has been developed by Dr. M Besora at the URV with the help of Dr. JI Mujika and Gonzalo D. Núñez.

## Citing NEST 

G. Zuccarello, L.J. Nannini, A. Arroyo-Bondía, N. Fincias, I. Arranz, A.H. Pérez-Jimeno, M. Peeters, I. Martín-Torres, A. Sadurní, V. García, Y. Wang, M.S. Kirillova, M. Montesinos-Magraner, U. Caniparoli, G.D. Núñez, F. Maseras, M. Besora, I. Escofet, A.M. Echavarren. Enantioselective Catalysis with Pyrrolidinyl Gold(I) Complexes: DFT and NEST Analysis of the Chiral Binding Pocket. To be published

