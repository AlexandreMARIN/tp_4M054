# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:25:11 2019

@author: 3200661
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri#pour les triangulations
from mpl_toolkits.mplot3d import Axes3D

def square_gmsh_file(R, lc, name):
    
    f = open(name+".geo", "w")
    f.write("R="+str(R)+";\nlc = "+str(lc)+";\n")
    dim1 = ["0", "R", "R", "0"]
    dim2 = ["0", "0", "R", "R"]
    for i in range(4):
        f.write("Point("+str(i+1)+") = { "+dim1[i]+" , "+dim2[i]+" , 0 , lc};\n")
    for i in range(4):
        f.write("Line("+str(i+1)+")= {"+str(i+1)+", "+str(((i+1)%4)+1)+"};\n")
    f.write("Line Loop(0) = {1, 2, 3, 4};\nPlane Surface(0) = {0};\nPhysical Line(0) = {0};")
    f.close()




def read_gmsh_file(file):
    """
    
    
    
    
    
    procédure:
        -sauter $MeshFormat
        -lire infos sur noeuds après $Nodes
        -lire infos sur elem après $Elements
        -attention aux tags + type éléments
            - 15 : point
            - 1 : ligne
            - 2 : triangle
    """
    f = open(file, "r")
    f.readline()
    f.readline()
    f.readline()
    
    f.readline()
    nodesnb = int(f.readline())
    nodes = np.zeros((nodesnb, 2))
    nodes_ind = np.array(nodesnb, dtype=np.int32)
    for i in range(nodesnb):
        aux = f.readline().split(" ")
        nodes[i][0] = float(aux[1])
        nodes[i][1] = float(aux[2])
        nodes_ind[i] = int(aux[0])

    f.readline()
    f.readline()
    elemnb = int(f.readline())#nombre d'arêtes sur le bord+nombre de triangles
    tab = np.empty((elemnb, 3), dtype=np.int32)
    for i in range(elemnb):
        aux = f.readline().split(" ")
        tab[i][0] = int(aux[1])
        tab[i][1] = int(aux[2])
        tab[i][2] = int(aux[3])
    f.close()
    return tab


#test

square_gmsh_file(1, 0.1, "square")
subprocess.run(["/Program Files/gmsh-4.1.0/gmsh", "square.geo", "-3"])