"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import pickle
import numpy as np
import scipy.sparse as sps

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

def extract_multilayer_degree(N,M,sparse_adjacency):
    """ Extract the multidegree sequence from the formatted sparse supra-adjacency matrix
            Input:
            N -- number of nodes
            M -- number of layers
            sparse_adjacency -- formatted sparse supra-adjacency matrix
            Output:
            l_Kx -- nodewise multidegree sequence
            l_Ky -- layerwise multidegree sequence
            """
    l_Kx = [0]*N
    l_Ky = [0]*M
    sparse_transpose = sparse_adjacency.transpose()
    for k in range(M):
        for i in range(N):
            degree_i_in_k = sparse_adjacency[i+(k*N)].nnz + sparse_transpose[i+(k*N)].nnz # degree of node i in layer space k
            l_Kx[i] += degree_i_in_k #add to node i k times
            l_Ky[k] += degree_i_in_k #add to node k i times
    return l_Kx,l_Ky


def formatting_supra_adjacency(N,M,sparse_matrix):
    """Format the  undirected supra-adjacency matrix in order to store the minimum information
    i.e. links information is stored in the upper diagonal of each block matrix
            Input:
            N -- number of nodes
            M -- number of layers
            sparse_matrix -- sparse supra-adjacency matrix
            Output:
            s -- formatted sparse supra-adjacency matrix
            """
    s = sps.lil_matrix((N * M, N * M), dtype=np.int8)
    l_nonzero = sparse_matrix.nonzero()
    for n in range(sparse_matrix.nnz):
        row,col = l_nonzero[0][n],l_nonzero[1][n]
        k,l,i,j = to_tensor_index(N,[row,col])
        if i < j:
            s[row,col] = 1
        if i == j and k<=l:
            s[row, col] = 1
    return s


if __name__ == "__main__":
    #LOAD DATA
    ##PARAMETERS SPACE
    load_path_1 = "Urban-multiplex-networks-main/adjacency_matrices/Duesseldorf_Aintra_unweighted.npz" #load intralayer matrices
    load_path_2 = "Urban-multiplex-networks-main/adjacency_matrices/Duesseldorf_Ainter_unweighted.npz" #load interlayer matrices
    save_path = "multilayer_network/german_transport/Duesseldorf_multiplex.dat" #where the data will be saved

    N = 1878 # number of nodes (needed to be known beforehand)
    M = 104 # number of layers (needed to be known beforehand)
    ##
    Aintra = sps.load_npz(load_path_1)
    Ainter = sps.load_npz(load_path_2)

    A = Ainter+Aintra
    A_formatted = formatting_supra_adjacency(N,M,A)

    lK_x,lK_y = extract_multilayer_degree(N,M,A_formatted)

    print("Total edges: ",sum(lK_x)/2)
    ##Save data
    dico = {}
    dico['sparse'] = A_formatted
    dico['N'] = N
    dico['M'] = M
    dico['l_K'] = lK_x
    dico['l_Ktranspose'] = lK_y
    # with open(save_path, 'wb') as f:
    #     pickle.dump(dico, f)



