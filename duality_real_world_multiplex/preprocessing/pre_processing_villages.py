"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import scipy.sparse as sps
import pandas as pd

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

def undirect_sparse_multiplex(N,M,sparse_matrix):
    """ Undirect a sparse adjacency matrix of a multiplex network
    Input:
    N -- number of nodes
    M -- number of layers
    sparse_matrix -- sparse supra-adjacency matrix
    Output:
    s -- undirected sparse supra-adjacency matrix"""
    s = sps.lil_matrix((N * M, N * M), dtype=np.int8)
    l_nonzero = sparse_matrix.nonzero()
    for n in range(sparse_matrix.nnz):
        row, col = l_nonzero[0][n], l_nonzero[1][n]
        k, l, i, j = to_tensor_index(N, [row, col])
        [row1,col1] = to_supra_index(N, [k,l,j,i])
        s[row, col] = 1
        s[row1, col1] = 1
    return s

if __name__ == "__main__":
    ##PARAMETERS SPACE
    load_path = "Uganda_villages/" #where the multiplex networks are
    save_path = "multilayer_network/Uganda_villages" #where to save them

    l_directory = [("friendship-"+str(i),"health-advice_"+str(i)) for i in range(1,18)]
    for i in range(len(l_directory)):
        df_friendship = pd.read_csv (load_path+l_directory[i][0]+"/nodes.csv")
        df_health = pd.read_csv (load_path+l_directory[i][1]+"/nodes.csv")
        l_index_to_friendship = []
        l_index_to_health = []
        name_friendship = list(df_friendship.loc[:," name"])
        name_health = list(df_health.loc[:," name"])

        names = list(set(name_friendship+name_health))

        for k in range(len(name_health)): #collect the indexed names of the nodes in health advice layer
            x = name_health[k]
            l_index_to_health.append(names.index(x) if x in names else None)

        for k in range(len(name_friendship)): #collect the indexed names of the nodes in friendship layer
            x = name_friendship[k]
            l_index_to_friendship.append(names.index(x) if x in names else None)

        N=len(names)
        M=2
        #for file in l_directory:

        # N = 246
        # M = 3
        ##
        file = load_path+l_directory[i][0]+"/edges.csv" #load data of the friendship layer
        A = sps.lil_matrix((N*M,N*M),dtype=np.int8)
        with open(file,'r') as f: #read the edge csv file
            while True:
                next_line = f.readline()
                if not next_line:
                    break
                target,source = next_line.split(",")
                nodeID_1 = int(target)
                nodeID_2 = int(source)
                nodeID_1 = l_index_to_friendship[nodeID_1] #associate index of node to the node itself
                nodeID_2 = l_index_to_friendship[nodeID_2] #associate index of node to the node itself
                row,col = to_supra_index(N,[0,0,nodeID_1,nodeID_2])
                A[row,col] = 1

        file = load_path+l_directory[i][1]+"/edges.csv" #load data of the health advice layer
        with open(file,'r') as f:
            while True:
                next_line = f.readline()
                if not next_line:
                    break
                target,source = next_line.split(",")
                nodeID_1 = int(target)
                nodeID_2 = int(source)
                nodeID_1 = l_index_to_health[nodeID_1]
                nodeID_2 = l_index_to_health[nodeID_2]
                row,col = to_supra_index(N,[1,1,nodeID_1,nodeID_2])
                A[row,col] = 1

        #
        A_undirected = undirect_sparse_multiplex(N,M,A)
        A_formatted = formatting_supra_adjacency(N,M,A_undirected)
        #

        lK_x,lK_y = extract_multilayer_degree(N,M,A_formatted)
        ## Save data
        dico = {}
        dico['sparse'] = A_formatted
        dico['N'] = N
        dico['M'] = M
        dico['l_K'] = lK_x
        dico['l_Ktranspose'] = lK_y
        # with open(save_path +'/village_'+str(i+1)+'.dat', 'wb') as f:
        #     pickle.dump(dico, f)