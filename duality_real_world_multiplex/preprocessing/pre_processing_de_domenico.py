"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import pickle
import numpy as np
import scipy.sparse as sps

def to_supra_index(N,uplet):
    """Convert tensorial indices (layer,layer,node,node) into flattened matrix indices (row,col)
        Input:
        N -- number of nodes
        uplet -- tensor indices (layer,layer,node,node)
        Output:
        row, col -- matrix indices (row,col)"""
    k, l, i, j = uplet[0],uplet[1],uplet[2],uplet[3]
    row = k*N +i
    col = l*N +j
    return row,col

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
    #Parameters
    load_path_1 = "EUAir_Multiplex_Transport/Dataset/EUAirTransportation_multiplex.edges" #where to get the edges
    save_path = "EUAir_multiplex" #where to save the data
    N = 450 #number of nodes (need to be known)
    M = 37 #number of layers
    ##RUN
    A = sps.lil_matrix((N*M,N*M),dtype=np.int8)
    with open(load_path_1,'r') as f:
        while True:
            next_line = f.readline()
            if not next_line:
                break
            layerID,nodeID_1,nodeID_2,weight = next_line.split(" ")
            layerID = int(layerID)-1 ## cast into integer and take out 1 to fit in sparse matrix from 0 to M-1
            nodeID_1 = int(nodeID_1) - 1
            nodeID_2 = int(nodeID_2) - 1
            row,col = to_supra_index(N,[layerID,layerID,nodeID_1,nodeID_2])
            A[row,col] = 1 #the supra-adjacency matrix is unweighted by this operation.

    A_formatted = formatting_supra_adjacency(N,M,A)
    #
    lK_x,lK_y = extract_multilayer_degree(N,M,A_formatted)
    #
    ## Save data
    dico = {}
    dico['sparse'] = A_formatted
    dico['N'] = N
    dico['M'] = M
    dico['l_K'] = lK_x
    dico['l_Ktranspose'] = lK_y
    # with open(save_path + '.dat', 'wb') as f:
    #     pickle.dump(dico, f)