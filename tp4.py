# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri#pour les triangulations
from mpl_toolkits.mplot3d import Axes3D


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
    Matrice de masse pour un triangle
    exo4 q2
    """
    M_T = np.array([[1.0/6, 1.0/24, 1.0/24], [1.0/24, 1.0/6, 1.0/24], [1.0/24, 1.0/24, 1.0/6]])
    area = 0.5*abs( (s1[0]-s2[0])*(s2[1]-s3[1]) - (s1[1]-s2[1])*(s2[0]-s3[0]) )

    M_T *= area
    return M_T


def RigElem(s1, s2, s3):
    """
    Matrice de rigidité pour un triangle
    """
    K_T = np.empty( (3, 3) )
    area = 0.5*abs( (s1[0]-s2[0])*(s2[1]-s3[1]) - (s1[1]-s2[1])*(s2[0]-s3[0]) )
    vert = np.array([np.array(s1), np.array(s2), np.array(s3)])
    for i in range(3):
        for j in range(3):
            K_T[i, j] = np.transpose((vert[(i+2)%3]-vert[(i+1)%3])).dot((vert[(j+2)%3]-vert[(j+1)%3]))
    K_T /= 4.0*area
    return K_T



def build_M(nodes, triangles):
    """
    Matrice de masse
    """
    M = np.zeros( (nodes.shape[0], nodes.shape[0]) )
    for q in range(triangles.shape[0]):
        aux = MassElem(nodes[triangles[q, 0]], nodes[triangles[q, 1]], nodes[triangles[q, 2]])
        for l in range(3):
            for m in range(3):
                M[triangles[q, m], triangles[q, l]] += aux[l, m]

    return M


def build_K(nodes, triangles):
    """
    Matrice de rigidité du système
    """
    K = np.zeros( (nodes.shape[0], nodes.shape[0]) )
    for q in range(triangles.shape[0]):
        aux = RigElem(nodes[triangles[q, 0]], nodes[triangles[q, 1]], nodes[triangles[q, 2]])
        for l in range(3):
            for m in range(3):
                K[triangles[q, m], triangles[q, l]] += aux[l, m]

    return K


def build_F(M, nodes, f):
    """
    construit le second membre
    M doit être la matrice de masse
    """
    vec = np.empty( (nodes.shape[0]) )
    for i in range(nodes.shape[0]):
        vec[i] = f(nodes[i][0], nodes[i][1])
    return M.dot(vec)

#tests

#demander à l'utilisateur le nombre de points en horizontal/vertical
hnb = int(input("Give an integer : "))
vnb = int(input("Give an integer : "))

create_rect_mesh("rect.msh", hnb, vnb, 1, 1)
#vérification du maillage
draw_mesh("rect.msh")

nodes = read_nodes_file("rect.msh")#alias (s_i)_i
print("nodes")
print(nodes)

print("triangles")
elements = read_elements_file("rect.msh")#alias I
print(elements)

#calcul de la matrice du système et du second membre
M = build_M(nodes, elements)
print("Mass")
print(M)

K = build_K(nodes, elements)
print("Rigidité :")
print(K)

#définition de deux fonctions u et f vérifiant l'équation (4.1)
def u(x, y):
    return np.cos(np.pi*x)*np.cos(np.pi*y)

def f(x, y):
    return (2*(np.pi**2)+1)*u(x, y)


F = build_F(M, nodes, f)
print("Second membre F :")
print(F)

print("A := M + K =")
print(M+K)

#calcul de la solution
U = np.linalg.solve(M+K, F)

print("U =")
print(U)

X = nodes[:, 0]
Y = nodes[:, 1]


#reconstruction de la solution à l'aide des triangulations
triang = mtri.Triangulation(X, Y, elements)
fig = plt.figure()
plt.title("$u(x,\, y)=\cos(\pi x)\,\cos(\pi y)$")
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.set_title("approximation")
ax.plot_trisurf(X, Y, U, triangles=elements, linewidth=0.2, antialiased=True)





ax = fig.add_subplot(1, 2, 2, projection='3d')
gridX = np.arange(0, 1, 0.05)
gridY = np.arange(0, 1, 0.05)
gridX, gridY = np.meshgrid(gridX, gridY)

U_exact = np.empty( gridX.shape )
for i in range(U_exact.shape[0]):
    for j in range(U_exact.shape[1]):
        U_exact[i, j] = u(gridX[i, j], gridY[i, j])
ax.plot_surface(gridX, gridY, U_exact, linewidth=0.2, antialiased=True)

ax.set_title("solution exacte")

plt.show()
