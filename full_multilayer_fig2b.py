import os
import pickle
import numpy as np
from math import fsum

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

def matrix_prepocessing(s_B): #
    """put all the values on top-right corner of the supra_adjacency matrix s_B
    Input:
    s_B -- supra-adjacency matrix
    Output:
    s_B -- framed supra-adjacency matrix"""
    row = s_B.nonzero()[0]
    col = s_B.nonzero()[1]
    A = s_B.toarray()
    for x in range(len(row)):
        a = row[x]
        b = col[x]
        k,l,i,j = to_tensor_index(N,[a,b])
        if k>l:
            A[b][a] = 1
            A[a][b] = 0
    return A

if __name__ == "__main__":
    # Load data
    print("Begin")
    name_path = "mock_data/full_N_200_M_200/"
    l_load_path = os.listdir(name_path)
    l_R = np.linspace(0,1,25)
    counter = 0
    l_proba_event = [[0.5, 0.5, 0], [0.25, 0.75, 0], [0.75, 0.25, 0]] #combination of proba to compute
    l_dD0 = np.zeros([len(l_proba_event),len(l_R)]) #list of distances in nodewise
    l_dDS = np.zeros([len(l_proba_event),len(l_R)]) #list of distances in layerwise
    l_density = np.zeros([len(l_proba_event),len(l_R)])
    for restore_path in l_load_path: #for each supra-adjacency matrix, compute its theoretical distances for a given r anf for all combination of rewiring event
        print(counter)
        with open(name_path+restore_path,'rb') as f:
            parameter = pickle.load(f) #we load the data file containing almost all the parameters already
        # parameters characterizing the random multilayer networks
        s_A = parameter['sparse']
        multi_degree = parameter['l_K']
        multi_degree_DS = parameter['l_Ktranspose']
        N = parameter['N']
        M = parameter['M']
       # l_density.append(parameter['density'])
        R = l_R[counter]*s_A.nnz
        counter += 1
        n = R
        A = matrix_prepocessing(s_A)
        number_edges = s_A.nnz
        C = np.zeros(shape=(N*M,N*M))
        l_K_C = [0]*N
        l_K_C_transpose = [0]*M
        number_pair = M * N * (M * N) / 2  # full multilayer case
        pair_links = M * M   # full multilayer case
        pair_layers = N * N / 2
        a = [0]*M
        l = [0]*M
        pair_tel = number_pair-pair_layers-pair_links+1
        density = number_edges/number_pair
        counter_proba = 0
        for p_space,p_edge,p_tel in l_proba_event:
            l_density[counter_proba][counter] = parameter['density']
            for k in range(M):
                for k in range(M):
                    l_K_C_transpose[k] = multi_degree_DS[k] + n * p_edge * (
                                2 * M / (pair_links - 1) - (pair_links * multi_degree_DS[k]) / (
                                    number_edges * (pair_links - 1)))
                    l_K_C_transpose[k] += n * p_tel * (
                            2 * M * (pair_layers - 1) / pair_tel - multi_degree_DS[k] * (number_pair - pair_links) / (
                            pair_tel * number_edges))
                for i in range(N):
                    l_K_C[i] = multi_degree[i] + n * p_space * (
                            (N) / (pair_layers - 1) - (pair_layers * multi_degree[i]) / (number_edges * (pair_layers - 1)))
                    # print(l_K_C[i])
                    l_K_C[i] += n * p_tel * (
                            (N) * (pair_links - 1) / pair_tel - multi_degree[i] * (number_pair - pair_layers) / (
                            pair_tel * number_edges))

            if fsum(multi_degree_DS) - fsum(l_K_C_transpose) > 0.1 or fsum(multi_degree) - fsum(l_K_C) > 0.1:
                print("Error for ", p_space,p_edge,p_tel)
                print(fsum(multi_degree_DS) - fsum(l_K_C_transpose))
            #print(l_K_C_transpose)
            delta = np.subtract(l_K_C, multi_degree)
            l_dD0[counter_proba][counter] = np.sqrt(np.dot(delta, delta))
            delta_DS = np.subtract(l_K_C_transpose, multi_degree_DS)
            l_dDS[counter_proba][counter] = np.sqrt(np.dot(delta_DS, delta_DS))
            counter_proba +=1

    dico = {}
    print(l_dD0,l_dDS)
    dico['sparse'] = s_A
    dico['l_distance_K0'] = l_dD0
    dico['l_distance_KDS'] = l_dDS
    dico['N'] = N
    dico['M'] = M
    dico['density'] = l_density
    dico['l_R'] = l_R
    dico['l_proba_event'] = l_proba_event
    save_path = "mock_data/figure_2_full_model_2_N_"+str(N)+"_M_"+str(M)+"_25.dat"
    with open(save_path, 'wb') as f:
        pickle.dump(dico, f)