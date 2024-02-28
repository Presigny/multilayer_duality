"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import numpy.random as rd
import pickle
from math import fsum
from time import perf_counter
import itertools

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
    return (row,col)

def edges_gen(N,s,layer_index,proba,replica=False):
    """Generate the random links in the layer_index[0],layer_index[1] according to the probability proba
        Input:
        N -- number of nodes
        s -- sparse supra-adjacency matrix of the being generated random sparse multilayer network
        layer_index --  (k,l) duplet of the ongoing layer index (k=l, intralayer, otherwise interlayer)
        proba -- probability to make a link in k,l (inter)layer
        replica -- if True, it adds  replica links at random according to the "proba" probabilities
        Output:
        l_K -- contribution of the k,l (inter)layer to the node multidegree sequence
        Ktranspose -- number of links to add to the k,l (inter)layer"""
    k,l = layer_index[0],layer_index[1]
    j = 0 #number of links to add to the k,l (inter)layer
    l_K = [0]*N
    edges = itertools.combinations(range(N),2) #all unique possible combination of links between nodes (iterable)
    for e in edges:
        if rd.random() < proba: #test the probability of appeareance
            j += 1
            row,col = to_supra_index(N,[k,l,e[0],e[1]])
            s[row,col] = 1 #place the link in the sparse supara-adjacency matrix
            l_K[e[0]] += 1 #add one to the multidegree of node e[0]
            l_K[e[1]] += 1
    if replica and k<=l: #add random replica links
        for i in range(N):
            if rd.random() < proba:
                j +=1
                row, col = to_supra_index(N,[k,l,i,i])
                s[row,col] = 1
                l_K[i] += 2
    K_transpose = j #number of links to add to the k,l (inter)layer
    return l_K,K_transpose

def sparse_generator(N,M,l_activity,replica=True,interlayer=True):
    """Generate the random sparse multilayer networks and its related node and layer multidegree sequence
        note 1: the generated sparse supra-adjacency matrix needs to be symmetrized to really represent a supra-adjacency matrix
        note 2: the replica links are not generated in this version,(see replica=False)
            Input:
            N -- number of nodes
            M -- number of layers
            l_activity --  activity (square root of probability to make link) of each layer
            replica -- boolean, if True presence of multiplex links (but randomly drawn, not full)
            interlayer -- boolean, if True presence of interlayer links
            Output:
            s -- sparse matrix (LIL format) of the generated network
            l_K -- multidegree sequence of the nodes in the generated network
            l_Ktranspose -- multidegree sequence of the layers in the generated network"""
    # step 0 (intralayer)
    l_Ktranspose = [0]*M
    l_K = [0]*N
    s = sps.lil_matrix((N*M,N*M),dtype=np.int8)
    if interlayer: #if it is a full multilayer network
        for k in range(M):
            print("Line ", k)
            for l in range(M):
                if k==l: #intermediate degree distributions
                    l_Kinter, Ktranspose_inter = edges_gen(N,s,[l,k], l_activity[k] * l_activity[k], replica=False)
                    l_Ktranspose[k] += 2*Ktranspose_inter # when within a layer a given link counts twice for the layer degree
                    l_K = [l_K[i] + l_Kinter[i] for i in range(len(l_K))]
                else:
                    l_Kinter, Ktranspose_inter = edges_gen(N,s,[l,k], l_activity[k] * l_activity[l], replica=replica)
                    l_Ktranspose[k] += Ktranspose_inter #the number of added links in itnerlayer goes to the k and l layer
                    l_Ktranspose[l] += Ktranspose_inter
                    l_K = [l_K[i] + l_Kinter[i] for i in range(len(l_K))]
    else: #if it is a multiplex network
        for k in range(M):
            print("Line ", k)
            l_Kinter, Ktranspose_inter = edges_gen(N,s,[k,k], l_activity[k] * l_activity[k], replica=False)
            l_Ktranspose[k] += 2*Ktranspose_inter # when within a layer a given link counts twice for the layer degree
            l_K = [l_K[i] + l_Kinter[i] for i in range(len(l_K))]
    return s,l_K,l_Ktranspose


if __name__ == "__main__":
    #Parameters
    n = 1 #number of multilayer network to generate
    N = 100 #number of desired nodes
    M = 5 #number of desired layers
    density = 0.05 #desired global density/probability parameter
    interlayer = True  # presence of interlayer links in the generated network
    multiplex_random = True  # presence of multiplex links (but randomly drawn, not full)
    multiplex = False  # presence of all multiplex links if True
    proba_set = [np.sqrt(density)] * M  #set the particular density of each layer
    for _ in range(n):
        path = "mock_data/mock_random_multilayer_network_"+str(n) #path where to save the generated multilayer network
        departure = perf_counter()
        s_A,l_K,l_Ktranspose = sparse_generator(N,M,proba_set,replica=multiplex_random,interlayer=interlayer)
        print("Time: ", perf_counter() - departure, " seconds")
        ## Save the data in dictionary
        dico = {}
        dico['sparse'] = s_A
        dico['N'] = N
        dico['M'] = M
        dico['l_K'] = l_K
        dico['density'] = density
        dico['l_Ktranspose'] = l_Ktranspose
        dico['inter'] = interlayer
        dico['mplex_random'] = multiplex_random
        # with open(path+'.dat','wb') as f:
        #     pickle.dump(dico,f)