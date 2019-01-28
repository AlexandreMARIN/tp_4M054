# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri#pour les triangulations



def read_lines(file):
    """
    lit et affiche les lignes d'un fichier
    """
    f = open(file, "r")
    for line in f:
        print(line, end='')
    f.close()


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
    tab = np.empty((elemnb, 3))
    for i in range(elemnb):
        aux = f.readline().split(" ")
        tab[i][0] = float(aux[1])
        tab[i][1] = float(aux[2])
        tab[i][2] = float(aux[3])
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


def draw_func_on_mesh(coords, triangles, col):
    """
    affiche les points de coords avec la couleur col qui est fonction de coords
    'triangles' est le tableau de triangles à considérer
    remarque : on peut renseigner une valeur dans col par triangle.
    """
    x = coords[:, 0]
    y = coords[:, 1]
    triang = mtri.Triangulation(x, y, triangles)
    plt.tripcolor(triang, col)
    plt.show()


#à retenir : for i, b in enumerate(ligns)

#tests


read_lines("maillage-tp1.msh")

print(read_nodes_file("maillage-tp1.msh"))

print(read_elements_file("maillage-tp1.msh"))


#création d'un maillage sur un rectangle puis affichage de celui-ci
create_rect_mesh("test.msh.txt", 30, 30, 50, 50)
coords = read_nodes_file("test.msh.txt")
triangles = read_elements_file("test.msh.txt")
x = coords[:, 0]
y = coords[:, 1]


triang = mtri.Triangulation(x, y, triangles)
plt.triplot(triang)

plt.show()


#ex 4 section 3

col = np.cos((x-25)*np.pi/25)*np.cos((y-25)*np.pi/25)

plt.tripcolor(triang, col)
plt.show()



