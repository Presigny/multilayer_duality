"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""

import numpy as np
import pickle
from random import choices

def to_supra_index(N,uplet):
    """Convert tensorial indices (layer,layer,node,node) into flattened matrix indices (row,col)
    Input:
    N -- number of nodes
    uplet -- tensor indices (layer,layer,node,node)
    Output:
    row, col -- matrix indices (row,col)"""
    k, l, i, j = uplet[0],uplet[1],uplet[2],uplet[3]
    row = k*N +i
    col =  l*N +j
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

def build_instructions(uplet,com):
    """Set which type of rewiring are allowed or not with a boolean instruction
    ex: com["be_interlayer] = False forbids the rewiring algorithm to rewire links in an interlayer space
        Input:
        uplet -- proposed tensor indices (layer,layer,node,node) to rewire in
        com -- dictionary with two keys ["be_interlayer"] and ["replica"]
        Output:
        boolean -- False if the proposed rewiring is forbidden """
    a,b,c,d = uplet[0], uplet[1], uplet[2], uplet[3]
    if com["be_interlayer"] == False: # CARE: present a matrix that has no interlayer links otherwise it runs during infinity
        if a!=b: #if we are in an interlayer space for the new_uplet...
            return False #... then forbids this change by False
    elif com["replica"] == False:
        if (a!=b and c == d): # we selected a replica as new_uplet
            return False # the change is forbidden
    return True

def select_intra(A,N,uplet,rng,com,trial):
    """Return the hole where 'uplet' is rewired according to the p_node event,
    which keeps the (inter)layer indices fixed
          Input:
          A -- supra-adjacency matrix
          N -- number of nodes
          uplet --  tensor indices (layer,layer,node,node) of the link to rewire
          rng --  random number generator
          com -- dictionary with two keys ["be_interlayer"] and ["replica"]
          trial -- number of proposed new rewired position
          Output:
          [k,l,a,b] -- Position where the link 'uplet' is rewired, it links node a in layer k to node b in layer l """
    k, l, i, j = uplet[0], uplet[1], uplet[2], uplet[3]  # the coordinates of the selected uplet
    rand_N = list(rng.integers(N,size=[trial,2]))
    for a,b in rand_N:
        avoid_self_loops = not((a == b) and (k == l)) # condition to avoid rewiring in self-loops
        conditions = build_instructions([k,l,a,b],com) and avoid_self_loops #and avoid_rewire_inter #and avoid_double_replica
        row,col = to_supra_index(N,[k,l,a,b])
        if A[row,col] == 0 and conditions:
            return [k, l, a, b]
    return None

def select_inter(A,N,M,uplet,rng, com,trial):
    """Return the hole where 'uplet' is rewired according to the p_layer event,
    which keeps the node indices fixed
          Input:
          A -- supra-adjacency matrix
          N -- number of nodes
          M -- number of layers
          uplet --  tensor indices (layer,layer,node,node) of the link to rewire
          rng --  random number generator
          com -- dictionary with two keys ["be_interlayer"] and ["replica"]
          trial -- number of proposed new rewired position
          Output:
          [a,b,i,j] -- Position where the link 'uplet' is rewired, it links node i in layer a to node j in layer b """
    k, l, i, j = uplet[0], uplet[1], uplet[2], uplet[3]  # the coordinates of the selected uplet
    rand_M = list(rng.integers(M, size=[trial, 2]))
    for a, b in rand_M:
        conditions = build_instructions([a, b, i, j], com)
        avoid_self_loops = not (a == b and i == j)  # condition to avoid rewiring in self-loops
        row, col = to_supra_index(N, [a, b, i, j])
        if A[row,col] == 0 and conditions and avoid_self_loops:
            return [a, b, i, j]
    return None

def select_tel(A,N,M,uplet,rng,com,trial):
    """Return the hole where 'uplet' is rewired according to the p_tel event,
    which changes at least one of the layer index and one one the node index
          Input:
          A -- supra-adjacency matrix
          N -- number of nodes
          M -- number of layers
          uplet --  tensor indices (layer,layer,node,node) of the link to rewire
          rng --  random number generator
          com -- dictionary with two keys ["be_interlayer"] and ["replica"]
          trial -- number of proposed new rewired position
          Output:
          [a,b,c,d] -- Position where the link 'uplet' is rewired, it links node c in layer a to node d in layer b """
    k, l, i, j = uplet[0], uplet[1], uplet[2], uplet[3]  # the coordinates of the selected uplet
    rand_M1 = list(rng.integers(M, size=[trial, 2]))
    rand_N1 = list(rng.integers(N, size=[trial, 2]))
    rand_NM = [[rand_M1[i][0],rand_M1[i][1],rand_N1[i][0],rand_N1[i][1]] for i in range(len(rand_M1))]
    for a,b,c,d in rand_NM:
        teleport_condition = ((a != k or b != l) and (c!=i or d!=j)) or ((a != l or b != k) and (c!=j or d!=i))  # minimal condition for the teleport event
        teleport_condition_2 = ((a != k or b != l) and (c!=j or d!=i)) or ((a != l or b != k) and (c!=i or d!=j)) # minimal condition for the teleport event
        if teleport_condition and teleport_condition_2:
            avoid_self_loops = not ((d == c) and (a == b))  # condition to avoid rewiring in self-loops
            avoid_rewire_inter = not ((c == i and d == j) or (c == j and d == i))  # condition to avoid making p_layer events
            row, col = to_supra_index(N, [a, b, c, d])
            conditions = build_instructions([a, b, c, d], com) and avoid_self_loops and avoid_rewire_inter
            if A[row,col] == 0 and conditions:
                return [a, b, c, d]
    return None


def stochastic_rewiring(s_A,multi_degree,multi_degree_DS,event_proba,R,dic_instructions,number_of_trial_per_event=1200):
    """Return the hole where 'uplet' is rewired according to the p_tel event,
    which changes at least one of the layer index and one the node index
          Input:
          s_A -- supra-adjacency matrix
          multi_degree -- node multidegree sequence
          multi_degree_DS -- layer multidegree sequence
          event_proba -- array [p_node,p_layer,p_tel] giving the probaiblity of each related rewiring event
          R -- proportion of links to rewire (from 0 to 1)
          dic_instructions -- dictionary with two keys ["be_interlayer"] and ["replica"]
          number_of_trial_per_event -- number of proposed new rewired position
          Output:
          s_B -- supra-adjacency matrix of the rewired multilayer network
          l_K -- node multidegree sequence of the rewired network
          l_Ktranspose -- layer multidegree sequence of the rewired network
          number_of_unrewire_edge -- number of links that are not rewired due to any problem """
    R = int(round((s_A.nnz * R))) # convert the proportion into number of links
    number_of_unrewire_edge = 0
    s_A = s_A.tocsr()
    s_A = s_A +s_A.T
    s_A = s_A.tolil()
    s_B = s_A.copy() #save the original matrix, s_B is the one that is rewired
    l_K = multi_degree.copy()
    l_Ktranspose = multi_degree_DS.copy()
    N = len(l_K) # number of nodes
    M=len(l_Ktranspose) #number of layers
    for r in range(R): #iteratively rewire R links
        row,col = s_B.nonzero() # row and col of all existing link in the network
        random_index = rng.integers(0, len(row), size=1) #select one link uniformly at random among all exisitng link
        while s_B[row[random_index][0],col[random_index][0]] == 2: # 2 indicates that the link were rewired already
            random_index = rng.integers(0,len(row), size=1)
        classic_uplet = [row[random_index][0], col[random_index][0]]  # always the first one of the random list
        uplet = to_tensor_index(N,classic_uplet)
        new_uplet = None
        security = 0
        while new_uplet is None:
            event = choices(["intra", "inter", "tel"], event_proba, k=1)[0] # randomly choose the rewiring event according to the probability in event_proba
            if event == "intra": #p_node rewiring event
                new_uplet = select_intra(s_B,N, uplet, rng, dic_instructions, number_of_trial_per_event)
            elif event == "inter": #p_inter rewiring event
                new_uplet = select_inter(s_B,N,M, uplet, rng, dic_instructions, number_of_trial_per_event)
            elif event == "tel": #p_tel rewiring event
                new_uplet = select_tel(s_B,N,M, uplet, rng, dic_instructions, number_of_trial_per_event)
            security += 1
            if security > 5:  # if the structure of the net makes it impossible to move the uplet so we keep it at its place
                new_uplet = uplet  # hoping it is not the case for a lot of nodes (normally not because sparse net)
                print("activate security")
                number_of_unrewire_edge += 1
        k, l, i, j = uplet[0], uplet[1], uplet[2], uplet[3]
        u, v, w, x = new_uplet[0], new_uplet[1], new_uplet[2], new_uplet[3]
        # Update the multidegree sequences with the rewiring k,l,i,j -> u,v,w,x
        l_Ktranspose[k] -= 1
        l_Ktranspose[l] -= 1
        l_K[i] -= 1
        l_K[j] -= 1
        l_Ktranspose[u] += 1
        l_Ktranspose[v] += 1
        l_K[x] += 1
        l_K[w] += 1
        classic_new_uplet = to_supra_index(N,new_uplet)
        # Update the supra-adjacency matrix with the rewiring k,l,i,j -> u,v,w,x
        s_B[classic_new_uplet[0],classic_new_uplet[1]] = 2 # 2 indicates to the function that this links was rewired already
        s_B[classic_new_uplet[1], classic_new_uplet[0]] = 2
        s_B[classic_uplet[0], classic_uplet[1]] = 0
        s_B[classic_uplet[1], classic_uplet[0]] = 0
    return s_B,l_K,l_Ktranspose,number_of_unrewire_edge

if __name__ == "__main__":
    # Parameters
    print("Begin")
    name_path = 'mock_data/mock_full_multilayer_N_20_M_20.dat'
    save_path = "//"
    with open(name_path,'rb') as f: #loading the data
        parameter = pickle.load(f)
    s_A = parameter['sparse'] #supra-adjacency amtrix of the network to be rewired
    multi_degree = parameter['l_K'] # node multidegree sequence
    multi_degree_DS = parameter['l_Ktranspose'] # layer multidegree sequence
    R = 1 # proportion of edges to rewire
    number_of_trial_per_event = 1200  # for each trial 12 draw are tested before new_iplet is declared None
    event_proba = [0, 0, 1] # [p_node,p_layer,p_tel] giving the probaiblity of each related rewiring event
    # idx=74
    # save_path_original = "mock_data/sparse_rewiring_data/Nodewise_rule_"+str(idx)+"_N_"+str(N)+"_M_"+str(M)+"_R_"+str(R)
    # save_path_DS = "mock_data/sparse_rewiring_data/Layerwise_rule_"+str(idx)+"_N_"+str(N)+"_M_"+str(M)+"_R_"+str(R)
    rng = np.random.default_rng()
    dic_instructions = {}
    dic_instructions["be_interlayer"] = True  # if False, we dont rewire into interlayer space, False is needed for multiplex networks)
    dic_instructions["replica"] = True  # if False, we can't rewire nodes in replic

## Running rewiring algorithm
    s_B,l_K,l_Ktranspose,number_unrewire = stochastic_rewiring(s_A,multi_degree,multi_degree_DS,
                                               event_proba,R,dic_instructions,number_of_trial_per_event=number_of_trial_per_event)
## Saving data
    dico = {}
    dico['sparse'] = s_B
    dico['l_K'] = l_K
    dico['l_Ktranspose'] = l_Ktranspose
    dico['idx'] = event_proba
    dico['N'] = len(l_K)
    dico['M'] = len(l_Ktranspose)
    dico['density'] = s_A.nnz
    dico['R'] = R
    dico['event_proba'] = event_proba
    dico['number_of_unrewired_edge'] = number_unrewire
    # with open(save_path+"EuAir_idx_"+str(idx), 'wb') as f:
    #     pickle.dump(dico, f)

## Saving sparse_rewiring_data

# with open(save_path_original, 'a') as f:
#     f.write(str(l_K)+"\n")
# with open(save_path_DS, 'a') as f:
#     f.write(str(l_Ktranspose)+"\n")