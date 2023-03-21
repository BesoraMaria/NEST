#!/usr/bin/env python
# coding: utf-8
# This program/app has been developed by Dr. M Besora at the URV with the help of Dr. JI Mujika and G.D. Nunez

import streamlit as st
    
import numpy as np
import math
from math import sqrt

st.title("NEST Volume Calculation")



#Setting a dictionary with atom radii

vdw = {
        "H":1.09, "He":1.40, "Li":1.82, "Be":1.53, "B":1.92, "C":1.70, "N":1.55, "O":1.52, "F":1.47, "Ne":1.54, "Na":2.27, "Mg":1.73, "Al":1.84, "Si":2.10, "P":1.80, "S":1.80, "Cl":1.75, "Ar":1.88, "K":2.75, "Ca":2.31, "Sc":2.15, "Ti":2.11, "V":2.07, "Cr":2.06, "Mn":2.05, "Fe":2.04, "Co":2.00, "Ni":1.63, "Cu":1.40, "Zn":1.39, "Ga":1.87, "Ge":2.11, "As":1.85, "Se":1.90, "Br":1.85, "Kr":2.02, "Rb":3.03, "Sr":2.49, "Y":2.32, "Zr":2.23, "Nb":2.18, "Mo":2.17, "Tc":2.16, "Ru":2.13, "Rh":2.10, "Pd":1.63, "Ag":1.72, "Cd":1.58, "In":1.93, "Sn":2.17, "Sb":2.06, "Te":2.06, "I":1.98, "Xe":2.16, "Cs":3.43, "Ba":2.68, "La":2.43, "Ce":2.42, "Pr":2.40, "Nd":2.39, "Pm":2.38, "Sm":2.36, "Eu":2.35, "Gd":2.34, "Tb":2.33, "Dy":2.31, "Ho":2.30, "Er":2.29, "Tm":2.27, "Yb":2.26, "Lu":2.24, "Hf":2.23, "Ta":2.22, "W":2.18, "Re":2.16, "Os":2.16, "Ir":2.13, "Pt":1.72, "Au":1.66, "Hg":1.55, "Tl":1.96, "Pb":2.02, "Bi":2.07, "Po":1.97, "At":2.02, "Rn":2.20, "Fr":3.48, "Ra":2.83, "Ac":2.47, "Th":2.45, "Pa":2.43, "U":1.86, "Np":2.39, "Pu":2.43, "Am":2.44, "Cm":2.45, "Bk":2.44, "Cf":2.45, "Es":2.45, "Fm":2.45, "Md":2.46, "No":2.46, "Lr":2.46,
}

#Defining calc_dist as a definition to measure distances from two 3D points
def calc_dist(c_x,c_y,c_z,a_x,a_y,a_z):
    d=sqrt((c_x-a_x)**2+(c_y-a_y)**2+(c_z-a_z)**2)
    return d

#Defining the angle between two vectors
def ang_2vec (A, B):
    
    ang_rad= np.arccos(np.dot(A,B)/(np.linalg.norm(A)*np.linalg.norm(B)))
    return(ang_rad)

#Defining empty variables
atoms=[]
coor_x=[]
coor_y=[]
coor_z=[]

#Uploading the xyz file to work with and reading all lines in

user_file = st.file_uploader("Choose a file")


if user_file is not None:
    all_lines = user_file.readlines()

with st.form(key="form1"):
    #Define atom numbers WITH IMPUT. Python starts counting from zero, but this will be arranged by the program few lines below.
    prim_atom = int(st.text_input("Number of origin atom (Metal):", 0))
    seg_atom = int(st.text_input("Number of atom in y axis (Central Ligand):", 0))
    ter_atom = int(st.text_input("Number of atom in z axis (out of xy plane):", 0))

    #Define the atoms to delete

    notinc = ([int(item) for item in st.text_input('Atoms excluded in V calculation (enter atom number):').split()] )

    #Define the size of the box with input

    size = float(st.text_input('NEST Size: 3,5 + length + 3,5. Give value for "length" (default=2):  ',2) or "2")

    #Multiplication factor on radii?
    factor = float(st.text_input('Scale the Bondi radii: 1 for no scale, typical scale 1.17, (default 1.17):',1.17) or "1.17")

    #Cell size, the smaller the more accurate but the more expensive
    cell_size = float(st.text_input('Grid for calculation of boxes (default 0.1):',0.1) or "0.1")


    run=st.form_submit_button("Run NEST with these data")

    if run:
        
        
        #deleting the two initial lines that in the xyz files do not correspond to atoms but atom numbers and title/empty
        del(all_lines[0])
        del(all_lines[0])
            
        #if there is a final empty line must be deleted
        #splitting the information of the coordinates in atoms, x y and z
        for line in all_lines:
            split=line.split() if line.strip() != ''
            atoms.append(split[0].decode())
            coor_x.append(float(split[1]))
            coor_y.append(float(split[2]))
            coor_z.append(float(split[3]))
        
        #setting the atom numbers to start counting from zero
        a= int(prim_atom) -1
        b= int(seg_atom) -1
        c= int(ter_atom) -1
        
        #defining some empty variables
        newcoor_x=[]
        newcoor_y=[]
        newcoor_z=[]
        radia=[]
        h1=[]
        ht=[]
            
        #centering coordinates to first atom (prim_atom, a)
        for i in range(len(coor_x)):
            newcoor_x.append(float(coor_x[i] - coor_x[a]))
            newcoor_y.append(float(coor_y[i] - coor_y[a]))
            newcoor_z.append(float(coor_z[i] - coor_z[a]))
        
        #setting second atom (seg_atom, b) in the y axis    
        pointy= calc_dist(newcoor_x[b],newcoor_y[b],newcoor_z[b],0,0,newcoor_z[b])
        
        p1 = (newcoor_x[b], newcoor_y[b], 0)
        p2 = (0, -pointy, 0)
        
        anglexy = ang_2vec(p1,p2)
        
        if ((int((newcoor_x[b]*np.cos(anglexy) - newcoor_y[b]*np.sin(anglexy))*1000) / 1000) != 0):
            anglexynew=2*(math.pi) - anglexy
            anglexy=anglexynew

        newcoor2_x = newcoor_x[:]
        newcoor2_y = newcoor_y[:]
        
        for i in range(len(coor_x)):
            newcoor2_x[i] = (newcoor_x[i]*np.cos(anglexy) - newcoor_y[i]*np.sin(anglexy))
            newcoor2_y[i] = (newcoor_x[i]*np.sin(anglexy) + newcoor_y[i]*np.cos(anglexy))
        
        #setting third atom(ter_atom,c) to define the z axis.
        pointy2= calc_dist(newcoor2_x[b],newcoor2_y[b],newcoor_z[b],0,0,0)
        p3=(newcoor2_x[b],newcoor2_y[b],newcoor_z[b])
        p4=(0,-pointy2,0)
        anglexyz = ang_2vec(p3,p4)

        if ((int((newcoor2_y[b]*np.sin(anglexyz) + newcoor_z[b]*np.cos(anglexyz))*1000) / 1000) != 0):
            anglexyznew=2*(math.pi) - anglexyz
            anglexyz=anglexyznew

        newcoor3_y = newcoor_y[:]
        newcoor2_z = newcoor_z[:]
        
        for i in range(len(coor_x)):
            newcoor3_y[i] = (newcoor2_y[i]*np.cos(anglexyz) - newcoor_z[i]*np.sin(anglexyz))
            newcoor2_z[i] = (newcoor2_y[i]*np.sin(anglexyz) + newcoor_z[i]*np.cos(anglexyz))
        
        pointy3= calc_dist(newcoor2_x[c],0,newcoor2_z[c],0,0,0)
        p5=(newcoor2_x[c],0,newcoor2_z[c])
        p6=(0,0,-pointy3)
        anglexz = ang_2vec(p5,p6)
        
        if ((int((newcoor2_x[c]*np.cos(anglexz) - newcoor2_z[c]*np.sin(anglexz)) * 1000) / 1000) != 0):
            anglexznew=2*(math.pi) - anglexz
            anglexz=anglexznew

        newcoor3_x = newcoor2_x[:]
        newcoor3_z = newcoor2_z[:]
        
        for i in range(len(coor_x)):
            newcoor3_x[i] = (newcoor2_x[i]*np.cos(anglexz) - newcoor2_z[i]*np.sin(anglexz))
            newcoor3_z[i] = (newcoor2_x[i]*np.sin(anglexz) + newcoor2_z[i]*np.cos(anglexz))   
        
        #Now the new alignement has been set.
        
        #Deleting atoms we don't want to take into account
        for i in range(len(notinc)):
            num=int(notinc[i])-1
            del(newcoor3_x[num])
            del(newcoor3_y[num])
            del(newcoor3_z[num])
            del(atoms[num])
        
        #creating two lists releated to radii, in the order of the atoms in the coordinates
        for n,name in enumerate(atoms):
            #creating an ordered list of radii in the order of the atoms in coordinates, to be used in the future.
            r_vdw=vdw[name]
            radia.append(factor*r_vdw)
            #creating a list of the maximum radii of ocupation to discern if the cell is more than half occupied or not.
            h1 = (radia[n]-math.sqrt(radia[n]**2 - (cell_size/2)**2)) 
            ht.append( h1 - (math.pi)*(h1**2)*(3*radia[n] - h1)/(3*(cell_size)**2))
        
        #creating a rectangularxsquared box of cells of given size, the size of x and z is set as 3.5 except for the y axis but it cand be changed to another number if needed or to a variable to be set by the user. Then it will be cutted to give a cilindrical shape with spherical ends
        x_min= -3.5
        x_max= 3.5
        y_min= -3.5
        y_max= 3.5 + size  
        z_min= -3.5
        z_max= 3.5
        
        #defining number of cells in each axis
        n_cells_x=int((x_max-x_min)/cell_size)
        n_cells_y=int((y_max-y_min)/cell_size)
        n_cells_z=int((z_max-z_min)/cell_size)
        
        #defining some empty variables
        n_z=0;n_y=0; n_x=0
        n=0
        free_cells=[]
        occup_cells=[]
        
        #setting a text to be opened in the screen and a box to keep a text changing until done.
        st.write('Take a little patience .. ' )
        screen_text = st.empty()
        
        #steping from a cell to another until done in all axis
        for n_z in range(n_cells_z):
            z_cell=z_min+(2*n_z+1)*(cell_size/2)
            #writing something on the screen to keep user confidence that there is some progress
            with screen_text.container():
                st.write('..', n_cells_z-n_z, '..')
            for n_y in range(n_cells_y):
                y_cell=y_min+(2*n_y+1)*(cell_size/2)
                for n_x in range(n_cells_x):
                    x_cell=x_min+(2*n_x+1)*(cell_size/2)
                    
                    #defining coordinates of the considered cell and defining "cell_occupied" it as "free" by default with the number "zero", number "1" is used for "occupied" and number "2" for cells out of the desired shape.
                    cell_coords=[x_cell, y_cell, z_cell]
                    cell_occupied=0
        
                    #computing the distance from the center of the cell to the y axis, the origin and the end of the cell. This is done in order to set the cilindrical shape and the spherical ends 
                    eixy = calc_dist(x_cell, y_cell, z_cell, 0, y_cell, 0)
                    eix0 = calc_dist(x_cell, y_cell, z_cell, 0, 0, 0)
                    eixend= calc_dist(x_cell, y_cell, z_cell, 0, size, 0)
                    #if the cell is out of the desired parameters / volume the "cell occupied" is set to number "2" that is given to cells out of the shape 
                    if ((y_cell <= 0.000 and eix0 >= 3.50001) or (y_cell > 0.000 and y_cell < size and eixy >=3.50001) or (y_cell >= size and eixend >= 3.50001)):
                        cell_occupied = 2
                    else:
                        #checking if the desired cell is occupied, if so number 1 is given.
                        n=0
                        for n in range(len(radia)):
                            radii=radia[n] + ht[n]
                            dist=calc_dist(x_cell, y_cell,z_cell,newcoor3_x[n],newcoor3_y[n],newcoor3_z[n])
                            if (radii > dist):
                                cell_occupied=1
                                break
                
                    if (cell_occupied==0):
                        free_cells.append(cell_coords)
                    elif (cell_occupied==1):
                        occup_cells.append(cell_coords)
                    n_x+=1
                n_y+=1
            n_z+=1
        
        #computing the total number of cells within the desired shape        
        real_total_cells=len(free_cells) + len(occup_cells)
        
        #computing the occupied volume
        NEST_volume= round((len(occup_cells))*100/real_total_cells,1)
        #printing results on screen
        st.write("NEST occupied volume:", NEST_volume,'%')
        
        a=0
        b=0
        
        
        #defining empty variables for quadrant division
        Q1_occ=[]
        Q2_occ=[]
        Q3_occ=[]
        Q4_occ=[]
        Q5_occ=[]
        Q1_free=[]
        Q2_free=[]
        Q3_free=[]
        Q4_free=[]
        Q5_free=[]
        
        #dividing cells into quadrants
        cell=0
        for cell in range(len(occup_cells)):
            if ((occup_cells[cell][0] >= 0) and (occup_cells[cell][1] >= 0) and (occup_cells[cell][2] >= 0)):
                Q1_occ.append(occup_cells[cell])
            elif ((occup_cells[cell][0] <= 0) and (occup_cells[cell][1] >= 0) and (occup_cells[cell][2] >= 0)):
                Q2_occ.append(occup_cells[cell])
            elif ((occup_cells[cell][0] >= 0) and (occup_cells[cell][1] >= 0) and (occup_cells[cell][2] <= 0)):
                Q3_occ.append(occup_cells[cell])
            elif ((occup_cells[cell][0] <= 0) and (occup_cells[cell][1] >= 0) and (occup_cells[cell][2] <= 0)):
                Q4_occ.append(occup_cells[cell])
            else: 
                Q5_occ.append(occup_cells[cell])
        
        cell=0
        
        for cell in range(len(free_cells)):
            if ((free_cells[cell][0] >= 0) and (free_cells[cell][1] >= 0) and (free_cells[cell][2] >= 0)):
                Q1_free.append(free_cells[cell])
            elif ((free_cells[cell][0] <= 0) and (free_cells[cell][1] >= 0) and (free_cells[cell][2] >= 0)):
                Q2_free.append(free_cells[cell])
            elif ((free_cells[cell][0] >= 0) and (free_cells[cell][1] >= 0) and (free_cells[cell][2] <= 0)):
                Q3_free.append(free_cells[cell])
            elif ((free_cells[cell][0] <= 0) and (free_cells[cell][1] >= 0) and (free_cells[cell][2] <= 0)):
                Q4_free.append(free_cells[cell])
            else:
                Q5_free.append(free_cells[cell])
        
        #computing quadrant results
        Q1occ=round(len(Q1_occ)/(len(Q1_occ)+len(Q1_free))*100,1)
        Q2occ=round(len(Q2_occ)/(len(Q2_occ)+len(Q2_free))*100,1)
        Q3occ=round(len(Q3_occ)/(len(Q3_occ)+len(Q3_free))*100,1)
        Q4occ=round(len(Q4_occ)/(len(Q4_occ)+len(Q4_free))*100,1)
        #printing results on screen
        st.write("Q1occ(x+,y+,z+)", Q1occ, "%")
        st.write("Q2occ(x-,y+,z+)", Q2occ, "%")
        st.write("Q3occ(x+,y+,z-)", Q3occ, "%")
        st.write("Q4occ(x-,y+,z-)", Q4occ, "%")

        
        
        #making outputfiles with geometries with all points and NEST
        print('Printing output geometry "geomNEST_'+str(user_file)+'" ready to open with jmol. Also four different files, named "geomNEST_Q?_'+str(user_file)+'", showing each quadrant.')
        out1 = open('geomNEST_'+user_file.name,'w')
        out1.write(str(len(newcoor3_x)+len(free_cells))+"\n")
        print('Total NEST occupancy:', NEST_volume, '%, Q1occ=', Q1occ,'%, Q2occ=', Q2occ,'%, Q3occ=', Q3occ,'%, Q4occ=', Q4occ,'%',file=out1)
        for a in range(len(newcoor3_x)):
            print(atoms[a], "%9.5f" % newcoor3_x[a], "%9.5f" % newcoor3_y[a], "%9.5f" % newcoor3_z[a],file=out1)
        for b in range(len(free_cells)):
            print('x', "%9.5f" % free_cells[b][0], "%9.5f" % free_cells[b][1], "%9.5f" % free_cells[b][2],file=out1)
        out1.close()
 
        
        a=0
        b=0
        out2= open('geomNEST_Q1_'+user_file.name,'w')
        out2.write(str(len(newcoor3_x)+len(Q1_free))+"\n")
        print('Results from NEST, Q1occ=', Q1occ, '%, Q2occ=', Q2occ, '%, Q3occ=', Q3occ,'%, Q4occ=', Q4occ,'%',file=out2)
        for a in range(len(newcoor3_x)):
            print(atoms[a], "%9.5f" % newcoor3_x[a], "%9.5f" % newcoor3_y[a], "%9.5f" % newcoor3_z[a],file=out2)
        for b in range(len(Q1_free)):
            print('x', "%9.5f" % Q1_free[b][0], "%9.5f" % Q1_free[b][1], "%9.5f" % Q1_free[b][2],file=out2)
        out2.close()
 
        
        a=0
        b=0
        out3= open('geomNEST_Q2_'+user_file.name,'w')
        out3.write(str(len(newcoor3_x)+len(Q2_free))+"\n")
        print('Results from NEST, Q1occ=',Q1occ,'%, Q2occ=',Q2occ,'%, Q3occ=', Q3occ,'%, Q4occ=', Q4occ,'%',file=out3)    
        for a in range(len(newcoor3_x)):
            print(atoms[a], "%9.5f" % newcoor3_x[a], "%9.5f" % newcoor3_y[a], "%9.5f" % newcoor3_z[a],file=out3)
        for b in range(len(Q2_free)):
            print('x', "%9.5f" % Q2_free[b][0], "%9.5f" % Q2_free[b][1], "%9.5f" % Q2_free[b][2],file=out3)
        out3.close()
        
        a=0
        b=0
        out4= open('geomNEST_Q3_'+user_file.name,'w')
        out4.write(str(len(newcoor3_x)+len(Q3_free))+"\n")
        print('Results from NEST, Q1occ=',Q1occ,'%, Q2occ=',Q2occ,'%, Q3occ=', Q3occ,'%, Q4occ=', Q4occ,'%',file=out4)
        for a in range(len(newcoor3_x)):
            print(atoms[a], "%9.5f" % newcoor3_x[a], "%9.5f" % newcoor3_y[a], "%9.5f" % newcoor3_z[a],file=out4)
        for b in range(len(Q3_free)):
            print('x', "%9.5f" % Q3_free[b][0], "%9.5f" % Q3_free[b][1], "%9.5f" % Q3_free[b][2],file=out4)
        out4.close()
        
        a=0
        b=0
        out5= open('geomNEST_Q4_'+user_file.name,'w')
        out5.write(str(len(newcoor3_x)+len(Q4_free))+"\n")
        print('Results from NEST, Q1occ=',Q1occ,'%, Q2occ=',Q2occ,'%, Q3occ=', Q3occ,'%, Q4occ=', Q4occ,'%',file=out5)
        for a in range(len(newcoor3_x)):
            print( atoms[a], "%9.5f" % newcoor3_x[a], "%9.5f" % newcoor3_y[a], "%9.5f" % newcoor3_z[a],file=out5)
        for b in range(len(Q4_free)):
            print('x', "%9.5f" % Q4_free[b][0], "%9.5f" % Q4_free[b][1], "%9.5f" % Q4_free[b][2],file=out5)
        out5.close()
        

XYZ= st.checkbox("Obtain xyz files of the NEST")
st.write("These files are very large and must be openned with a visualization program that allows to work with large amount of data. We recomment using Jmol. Jmol is available in Windows, Linux and Mac.")

if XYZ:
    with open('geomNEST_'+user_file.name,'r') as file:
        btn = st.download_button(
            label="Download NEST geometry",
            data=file,
            file_name='geomNEST_'+user_file.name
            )

    with open('geomNEST_Q1_'+user_file.name,'r') as file:
        btn = st.download_button(
            label="Download Q1",
            data=file,
            file_name='geomNEST_Q1_'+user_file.name
            )

    with open('geomNEST_Q2_'+user_file.name,'r') as file:
        btn = st.download_button(
            label="Download Q2",
            data=file,
            file_name='geomNEST_Q2_'+user_file.name
            )

    with open('geomNEST_Q3_'+user_file.name,'r') as file:
        btn = st.download_button(
            label="Download Q3",
            data=file,
            file_name='geomNEST_Q3_'+user_file.name
            )

    with open('geomNEST_Q4_'+user_file.name,'r') as file:
        btn = st.download_button(
            label="Download Q4",
            data=file,
            file_name='geomNEST_Q4_'+user_file.name
            )



