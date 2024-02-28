"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import os

def to_tensor_index(N,uplet):
    """Convert flattened matrix indices (row,col) into tensorial indices (layer,layer,node,node)
        Input:
        N -- number of nodes
        uplet -- matrix indices (row,col)
        Output:
        row, col -- tensor indices (layer,layer,node,node) """
    row,col = uplet[0],uplet[1]
    k = row // N
    l = col // N
    i = row % N
    j = col % N
    return k,l,i,j

def to_supra_index(N,uplet):
    """Convert tensorial indices (layer,layer,node,node) into flattened matrix indices (row,col)
            Input:
            N -- number of nodes
            uplet -- matrix indices (row,col)
            Output:
            row, col -- matrix indices (layer,layer,node,node) """
    k, l, i, j = uplet[0],uplet[1],uplet[2],uplet[3]
    row = k*N +i
    col =  l*N +j
    return row,col

## Load the path where the matrices to be reduces are
load_path = "database/sym_individual_matrix/"
N = 70
M = 77
## Select the layers to supress in the adjacency matrices

l_sup_layer = [0,12,25,38,51,63,76] #70 layers
#l_sup_layer = [ 0,  5, 10, 15 , 20,25, 30, 35, 41, 46,51, 56, 61, 66, 71,76] #61 layers
#l_sup_layer = np.linspace(0,76,23,dtype=int) #54 layers
#l_sup_layer = np.linspace(0,76,31,dtype=int) #46 layers
#l_sup_layer = np.linspace(0,76,39,dtype=int) #38 layers
#l_sup_layer = np.linspace(0,76,46,dtype=int) #31 layers
#l_sup_layer = np.linspace(0,76,54,dtype=int) #23 layers
#l_sup_layer = np.linspace(0,76,58,dtype=int) #19 layers
#l_sup_layer = np.linspace(0,76,62,dtype=int) #15 layers
#l_sup_layer = np.linspace(0,76,70,dtype=int) #7 layers
#l_sup_layer = np.linspace(0,76,74,dtype=int) #3 layers
number_layers = M-len(l_sup_layer)
# Select the path where to store the reduced matrix
dir_path = "database_reduced/%s_sym_individual_matrix/"%(number_layers)
os.mkdir(dir_path)

l_node = np.arange(N)
l_layer = np.arange(M)
l_load_path = os.listdir(load_path)

# Supression of the links in the matrices
for matrix_path in l_load_path:
    A = np.load(load_path + matrix_path)
    for k in l_layer:
        for l in l_sup_layer:
            for i in l_node:
                for j in l_node:
                    row,col = to_supra_index(N,[k,l,i,j])
                    A[row,col] = 0
                    A[col,row] = 0
    #Save data
    with open(dir_path+matrix_path, 'wb') as f:
        np.save(f,A)
    print(matrix_path, " DONE")