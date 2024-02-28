"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import os
import pickle

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


def multistrength(A,N):
    """ Compute the nodewise multistrength of supra-adjacency matrix A
        Input:
        A -- supra-adjacency matrix
        N -- number of nodes
        Output:
        multidegree -- nodewise multistrength of supra-adjacency matrix A
    """
    multistrength = [0]*N
    for row in range(len(A)):
        for col in range(row+1):
            k,l,i,j = to_tensor_index(N,[row,col]) # from i to j => out
            multistrength[i] += A[row,col]
            multistrength[j] += A[row,col]
    return multistrength

def multistrength_DS(A,N,M):
    """ Compute the layerwise multistrength of supra-adjacency matrix A
        Input:
        A -- supra-adjacency matrix
        N -- number of nodes
        M -- number of layers
        Output:
        multidegree -- layerwise multistrength of supra-adjacency matrix A
    """
    multistrength = [0] * M
    for row in range(len(A)):
        for col in range(row+1):
            k, l, i, j = to_tensor_index(N, [row, col])  # from i to j => out
            multistrength[k] += A[row, col]
            multistrength[l] += A[row, col]
    return multistrength

def group_by_type(load_path):
    """Separe the filenames of Patients and Subjects
        Input:
        load_path -- global path where patients and subjects names are mixed
        Output:
        l_group_suj -- subjects filename
        l_group_pat -- patients filename"""
    l_load_path = os.listdir(load_path)
    l_group_suj = []
    l_group_pat = []
    for name_1 in l_load_path:
        if name_1.count("PP") != 0:#detect patients filename
            l_group_pat.append(name_1)
        else:
            l_group_suj.append(name_1)
    return l_group_suj,l_group_pat

##
#load_path = "database/sym_individual_matrix/"
#load_path = "database/sym_avgSUJ_matrix/"
#dir_path = "database/sym_multistrength_SUJ/"
#dir_path = "database/sym_multistrength_avgSUJ/"

load_path = "database_reduced/70_sym_average_SUJ_matrix/"
dir_path = "database_reduced/70_sym_multistrength_avgSUJ/"
#os.mkdir(dir_path)
N = 70 #number of nodes

l_load_suj,l_load_pat = group_by_type(load_path)

#l_load_path = l_load_pat #select here only the patients for their multistrength to be computed
l_load_path = l_load_suj #select here only the subjects for their multistrength to be computed

for matrix_path in l_load_path:
    A = np.load(load_path+matrix_path)
    M = int(len(A)/N)
    k = multistrength(A,N)
    k_DS = multistrength_DS(A,N,M)
    dico = {}
    dico['multistrength'] = k
    dico['multistrength_DS'] = k_DS
    # save the nodewise/layerwise multistrength
    with open(dir_path+matrix_path+"_multistrength", 'wb') as f:
        pickle.dump(dico, f)
    print(matrix_path, " DONE")
