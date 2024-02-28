"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import pickle
import numpy as np
from math import fsum
import os

def lin_proba_set(res_tel,res_other,rounding = 3):
    """generate linearly a list of probability combinations that follow the rule: p_intra +p_inter <= 1 and p_tel = 1-(p_intra+p_inter)
            Input:
            res_tel -- approximate resolution of the p_tel parameter
            res_other -- approximate resolution of the p_node and p_layer parameters
            rounding -- generated probabilities are rounded at the ''rounding'' decimal
            Output:
            l_proba -- list of probability combination (numer of elements is approx.
            """
    l_proba = []
    step_tel = int(1/res_tel)
    step_other = int(1/res_other)
    ptel = 0
    pinter = 0
    for i in range(int(step_tel)):
        for j in range(int(step_other)):
            if pinter <= 1-ptel: #condition for those values to be probabilities ot select a given rewiring event.
                pintra = 1-ptel-pinter
                l_proba.append([round(pintra,rounding),round(pinter,rounding),round(ptel,rounding)])
            pinter += res_other
        pinter = 0
        ptel += res_tel
    ptel = 0 #this colum to add the variables p_intra=0
    for i in range(int(step_tel)):
            l_proba.append([0,round(1-ptel,rounding),round(ptel,rounding)])
            ptel += res_tel
    l_proba.append([0,0,1])
    return l_proba

def list_of_average_theoretical_distance_array_full(load_path,r,l_proba_event):
    """Compute the average theoretical multidegree distances for the input multilayer networks (cf. paper) in the nodewise and layerwise side
    for the input combinations of probabilities [p_node,p_layer,p_tel]. Note that each network in the input directory
    is associated to a given R proportion to rewire links (in list l_R).
              Input:
              load_path -- path to the directory where the full network is
              r -- proportion of links rewired in the network during the rewiring process
              l_proba_event -- list of probabilities of the node,layer and tel event triples (format - [p_node,p_layer,p_tel])
              Output:
              l_dD0 -- For each combination of probabilities, the list contains the theoretical distance in nodewise side between
                    the input multiplex networks and their theoretical averaged counterpart, format -> (len(l_proba_event),len(l_load_path))
              l_dDS -- For each combination of probabilities, the list contains the theoretical distance in layerwise side between
                    the input multiplex networks and their theoretical averaged counterpart, format ->  (len(l_proba_event),len(l_load_path))
              N -- number of nodes of networks in the directory
              M -- number of layers of networks in the directory (
              """
    l_load_path = os.listdir(load_path)
    l_R = np.linspace(0, 1, len(l_load_path))  # produce the list of R
    counter = 0  # index of the multiplex networks in distances list
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
        counter_proba = 0  # index of a given combination of proability in lsit l_proba_event
        number_edges = s_A.nnz
        l_K_C = [0] * N
        l_K_C_transpose = [0] * M
        number_pair = M * N * (M * N) / 2  # number of pairs full multilayer case
        pair_links = M * M  # number of possible links between a given pair of links full multilayer case
        pair_layers = N * N / 2 # number of possible links in a layer
        pair_tel = number_pair - pair_layers - pair_links + 1
        for p_space, p_edge, p_tel in l_proba_event: #loop over the combinations of probability in l_proba_event
            for k in range(M): # compute the theoretical average layer multidegree after rewiring
                l_K_C_transpose[k] = multi_degree_DS[k] + n * p_edge * (
                        2 * M / (pair_links - 1) - (pair_links * multi_degree_DS[k]) / (
                        number_edges * (pair_links - 1)))
                l_K_C_transpose[k] += n * p_tel * (
                        2 * M * (pair_layers - 1) / pair_tel - multi_degree_DS[k] * (number_pair - pair_links) / (
                        pair_tel * number_edges))
            for i in range(N): # compute the theoretical average nodewise multidegree after rewiring
                l_K_C[i] = multi_degree[i] + n * p_space * (
                        (N) / (pair_layers - 1) - (pair_layers * multi_degree[i]) / (number_edges * (pair_layers - 1)))
                l_K_C[i] += n * p_tel * (
                        (N) * (pair_links - 1) / pair_tel - multi_degree[i] * (number_pair - pair_layers) / (
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
            counter_proba += 1  # update the index of multiplex network
        counter += 1  # update the index of input combination probabilities
    return l_dD0,l_dDS,N,M


if __name__ == "__main__":
    ##Load the initial distance
    # Load data
    print("Begin")
    #Parameters
    res_other = 0.1 # approximate resolution of the p_node and p_layer parameters
    res_tel = 0.1 # approximate resolution of the p_tel parameter

    r=0.5 # proportion of links to rewire
    l_proba_event = lin_proba_set(res_tel, res_other)
    # (1/res_other + 1/res_tel)/2 + 1/p_tel gives approximatively the number of probability combination generated by lin_proba_set
    restore_path = "mock_data/full_N_200_M_200/"

    l_dD0,l_dDS,N,M= list_of_average_theoretical_distance_array_full(restore_path, r, l_proba_event)
    ##Save the data
    dico = {}
    dico['l_distance_K0'] = l_dD0
    dico['l_distance_KDS'] = l_dDS
    dico['N'] = N
    dico['M'] = M
    dico['R'] = r
    dico['l_proba_event'] = l_proba_event
    #save_path = "data/theoretical/a_full_model_2_N_"+str(N)+"_M_"+str(M)+"_R_"+str(r)+".dat"
    # with open(save_path, 'wb') as f:
    #     pickle.dump(dico, f)