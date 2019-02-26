# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri#pour les triangulations



def read_nodes_file(file):
    """
    lit les points constitutifs des triangles du maillage du fichier file
    """
    f = open(file, "r")
    f.readline()
    nodesnb = int(f.readline())
    tab = np.zeros((nodesnb, 3))
    for i in range(nodesnb):
        aux = f.readline().split(" ")
        tab[i][0] = float(aux[1])
        tab[i][1] = float(aux[2])
        if(len(aux)==4):
            tab[i][2] = float(aux[3])
    f.close()
    return tab

def read_elements_file(file):
    """
    lit les éléments d'un fichier de maillage du fichier file constituant les triangles
    """
    f = open(file, "r")
    f.readline()
    nodesnb = int(f.readline())
    for i in range(nodesnb):
        f.readline()
    f.readline()
    f.readline()
    elemnb = int(f.readline())
    tab = np.empty((elemnb, 3), dtype=np.int32)
    for i in range(elemnb):
        aux = f.readline().split(" ")
        tab[i][0] = int(aux[1])
        tab[i][1] = int(aux[2])
        tab[i][2] = int(aux[3])
    f.close()
    return tab



def create_rect_mesh(filename, hornb, vertnb, length, width):
    """
    ex 2 tp1
    crée le fichier filename contenant le maillage du rectangle
    de coin gauche inférieur l point (0, 0)
    de côtés length et width
    avec :
        hornb points selon l'horizontal
        vertnb points selon la verticale
    """

    f = open(filename, "w")
    f.write("$Noeuds\n")
    f.write(str(hornb*vertnb))
    f.write("\n")
    
    #noeuds
    hedge = float(length)/(hornb-1)
    vedge = float(width)/(vertnb-1)
    for i in range(vertnb):
        for j in range(hornb):
            f.write(str(i*hornb+j))
            f.write(" ")
            f.write(str(hedge*j))
            f.write(" ")
            f.write(str(vedge*i))
            f.write("\n")
    f.write("$FinNoeuds\n$Elements\n")
    trinb = (hornb-1)*(vertnb-1)*2
    f.write(str(trinb))
    f.write("\n")
    trinb //= 2#for calculating indices
    #éléments
    for i in range(vertnb-1):
        for j in range(hornb-1):
            f.write(str(i*hornb+j))
            f.write(" ")
            f.write(str(i*hornb+j))
            f.write(" ")
            f.write(str(i*hornb+j+1))
            f.write(" ")
            f.write(str((i+1)*hornb+j))
            f.write("\n")
            #second triangle
            f.write(str(trinb + i*hornb+j))
            f.write(" ")
            f.write(str((i+1)*hornb+j+1))
            f.write(" ")
            f.write(str(i*hornb+j+1))
            f.write(" ")
            f.write(str((i+1)*hornb+j))
            f.write("\n")
    f.write("$FinElements")
    f.close()


def draw_mesh(filename):
    """
    dessine le maillage décrit dans le fichier de nom 'filename'
    """
    coords = read_nodes_file(filename)
    triangles = read_elements_file(filename)
    x = coords[:, 0]
    y = coords[:, 1]
    
    triang = mtri.Triangulation(x, y, triangles)
    plt.triplot(triang)
    plt.show()



def MassElem(s1, s2, s3):
    """
    exo4 q2
    """
    M_T = np.array([[1.0/6, 1.0/24, 1.0/24], [1.0/24, 1.0/6, 1.0/24], [1.0/24, 1.0/24, 1.0/6]])
    area = 0.5*abs( (s1[0]-s2[0])*(s2[1]-s3[1]) - (s1[1]-s2[1])*(s2[0]-s3[0]) )

    M_T *= area
    return M_T


def RigElem(s1, s2, s3):
    pass




def build_M(nodes, triangles):
    M = np.zeros( (nodes.shape[0], nodes.shape[0]) )
    for q in range(triangles.shape[0]):
        aux = MassElem(nodes[triangles[q, 0]], nodes[triangles[q, 1]], nodes[triangles[q, 2]])
        for l in range(2):
            for m in range(2):
                M[triangles[q, m], triangles[q, l]] += aux[l, m]

    return M


def build_K():
    pass


#tests


create_rect_mesh("rect.msh", 2, 2, 1, 1)
draw_mesh("rect.msh")

nodes = read_nodes_file("rect.msh")#alias (s_i)_i
print("nodes")
print(nodes)

print("triangles")
elements = read_elements_file("rect.msh")#alias I
print(elements)


M = build_M(nodes, elements)
print("Mass")
print(M)
