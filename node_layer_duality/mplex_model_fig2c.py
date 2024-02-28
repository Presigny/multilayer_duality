"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import os
import pickle
import numpy as np
from math import fsum

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

def list_of_average_theoretical_distance_mplex(load_path,l_proba_event):
    """Compute the average theoretical multidegree distances for the input multiplex networks (cf. paper) in the nodewise and layerwise side
    for the input combinations of probabilities [p_node,p_layer,p_tel]. Note that each network in the input directory
    is associated to a given R proportion to rewire links (in list l_R).
              Input:
              load_path -- path to the directory where the multiplex networks are
              l_proba_event -- list of probabilities of the node,layer and tel event triples (format - [p_node,p_layer,p_tel])
              Output:
              l_dD0 -- For each combination of probabilities, the list contains the theoretical distance in nodewise side between
                    the input multiplex networks and their theoretical averaged counterpart, format -> (len(l_proba_event),len(l_load_path))
              l_dDS -- For each combination of probabilities, the list contains the theoretical distance in layerwise side between
                    the input multiplex networks and their theoretical averaged counterpart, format ->  (len(l_proba_event),len(l_load_path))
              N -- number of nodes of networks in the directory (the same for all networks)
              M -- number of layers of networks in the directory (the same for all networks)
              l_density -- list of density of the multiplex networks
              l_R -- list of proportions of links rewired in the network during the rewiring process"""
    l_load_path = os.listdir(load_path)
    l_R = np.linspace(0, 1, len(l_load_path)) #produce the list of R
    counter = 0 # index of the multiplex networks in distances list
    l_dD0 = np.zeros([len(l_proba_event), len(l_R)])
    l_dDS = np.zeros([len(l_proba_event), len(l_R)])
    l_density = []
    for restore_path in l_load_path:
        print(counter)
        with open(load_path + restore_path, 'rb') as f:
            parameter = pickle.load(f)  # we load the data file containing almost all the parameters already
        # parameters characterizing the random multilayer networks
        s_A = parameter['sparse']
        multi_degree = parameter['l_K']
        multi_degree_DS = parameter['l_Ktranspose']
        N = parameter['N']
        M = parameter['M']
        l_density.append(parameter['density'])
        R = l_R[counter] * s_A.nnz
        n = R
        number_edges = s_A.nnz
        l_K_C = [0] * N
        l_K_C_transpose = [0] * M
        number_pair = M * N * (N - 1) / 2  # number of pairs in the mplex case
        pair_links = M  # number of possible links between a given pair of links in mplex case
        pair_layers = N * (N - 1) / 2 # number of possible links in a layer
        pair_tel = number_pair - pair_layers - pair_links + 1
        counter_proba = 0 # index of a given combination of proability in lsit l_proba_event
        for p_space, p_edge, p_tel in l_proba_event: # loop over each input combination of probabilities
            for k in range(M): # compute the theoretical average layer multidegree after rewiring
                l_K_C_transpose[k] = multi_degree_DS[k] + n * p_edge * (
                        2 / (pair_links - 1) - (pair_links * multi_degree_DS[k]) / (number_edges * (pair_links - 1)))
                l_K_C_transpose[k] += n * p_tel * (
                        2 * (pair_layers - 1) / pair_tel - multi_degree_DS[k] * (number_pair - pair_links) / (
                        pair_tel * number_edges))
            for i in range(N): # compute the theoretical average nodewise multidegree after rewiring
                l_K_C[i] = multi_degree[i] + n * p_space * (
                        (N - 1) / (pair_layers - 1) - (pair_layers * multi_degree[i]) / (
                        number_edges * (pair_layers - 1)))
                l_K_C[i] += n * p_tel * (
                        (N - 1) * (pair_links - 1) / pair_tel - multi_degree[i] * (number_pair - pair_layers) / (
                        pair_tel * number_edges))
            if fsum(multi_degree_DS) - fsum(l_K_C_transpose) > 0.1 or fsum(multi_degree) - fsum(l_K_C) > 0.1:
                # check if the degree is conserved after rewiring
                print("Error for ", p_space, p_edge, p_tel)
                print(fsum(multi_degree_DS) - fsum(l_K_C_transpose))
            # compute the euclidean distance with node- and layer-wise multidegree
            delta = np.subtract(l_K_C, multi_degree)
            l_dD0[counter_proba][counter] = np.sqrt(np.dot(delta, delta))
            delta_DS = np.subtract(l_K_C_transpose, multi_degree_DS)
            l_dDS[counter_proba][counter] = np.sqrt(np.dot(delta_DS, delta_DS))
            counter_proba += 1 #update the index of multiplex network
        counter += 1 #update the index of input combination probabilities
    return l_dD0,l_dDS,N,M,l_density,l_R

if __name__ == "__main__":
    ##Load the initial distance
    # Load data
    print("Begin")
    # Parameters
    path = "mock_data/mplex_N_200_M_200/"  #path where to find the multiplex networks
    l_proba_event = [[0.5, 0.5, 0], [0.25, 0.75, 0], [0.75, 0.25, 0]] # list of combination of probabilities
    l_dD0,l_dDS,N,M,l_density,l_R = list_of_average_theoretical_distance_mplex(path, l_proba_event)

## Save data in a dictionary
    dico = {}
    dico['l_distance_K0'] = l_dD0
    dico['l_distance_KDS'] = l_dDS
    dico['N'] = N
    dico['M'] = M
    dico['density'] = l_density
    dico['l_R'] = l_R
    dico['l_proba_event'] = l_proba_event
    #Save parameters
    #save_path = "data/theoretical/figure_2_mplex_model_2_N_" + str(N) + "_M_" + str(M) + "_stat.dat" #path where to save the results of the script
    # with open(save_path, 'wb') as f:
    #     pickle.dump(dico, f)
    #print("Time spent: ", perf_counter() - departure, " seconds")