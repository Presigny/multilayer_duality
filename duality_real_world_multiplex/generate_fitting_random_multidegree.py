"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import pickle
import numpy as np

def distance(list1,list2):
   delta = np.subtract(list1, list2)
   return np.sqrt(np.dot(delta, delta))

def compute_theoretical_distance(multi_degree,multi_degree_DS,R,l_proba):
    """Compute the theoretical between nodewise/layerwise multidegree sequence of a multiplex network and
    its averaged counterpart as obtained by the stochastic rewiring algorithm
            Input:
            multi_degree -- nodewise multidegree sequence of a multiplex network
            multi_degree_DS -- layerwise multidegree sequence of a multiplex network
            R -- proportion of links to rewire
            l_proba -- combination of probability to use for the theoreticla stochastic rewiring
            Output:
            d_x --  theoretical  nodewise distance between the input multiplex network's multidegree sequence and its theoretical averaged counterpart
            d_y -- theoretical  layerwise distance between the input multiplex network's multidegree sequence and its theoretical averaged counterpart"""
    N = len(multi_degree)
    M = len(multi_degree_DS)
    R = R * np.sum(multi_degree_DS)/2
    p_space,p_edge,p_tel = l_proba[0],l_proba[1],l_proba[2]
    n = R
    number_edges = np.sum(multi_degree_DS)/2
    l_K_C = [0] * N
    l_K_C_transpose = [0] * M
    number_pair = M * N * (N - 1) / 2  # number of pairs in the mplex case
    pair_links = M  # number of possible links between a given pair of links in mplex case
    pair_layers = N * (N - 1) / 2  # number of possible links in a layer

    pair_tel = number_pair - pair_layers - pair_links + 1
    for k in range(M): # compute the theoretical average layer multidegree after rewiring
        l_K_C_transpose[k] = multi_degree_DS[k] + n * p_edge * (
                    2 / (pair_links - 1) - (pair_links * multi_degree_DS[k]) / (number_edges * (pair_links - 1)))
        l_K_C_transpose[k] += n * p_tel * (
                    2 * (pair_layers - 1) / pair_tel - multi_degree_DS[k] * (number_pair - pair_links) / (
                        pair_tel * number_edges))
    for i in range(N): # compute the theoretical average node multidegree after rewiring
        l_K_C[i] = multi_degree[i] + n * p_space * (
                    (N - 1) / (pair_layers - 1) - (pair_layers * multi_degree[i]) / (number_edges * (pair_layers - 1)))
        l_K_C[i] += n * p_tel * (
                    (N - 1) * (pair_links - 1) / pair_tel - multi_degree[i] * (number_pair - pair_layers) / (
                        pair_tel * number_edges))
    # compute the euclidean distance with node- and layer-wise multidegree
    d_x = distance(l_K_C, multi_degree)
    d_y = distance(l_K_C_transpose, multi_degree_DS)
    return d_x,d_y

def equalizer(l_new,l_ref):
    """Correct a randomly generated multidegree sequence to fit the total degree to the one of a reference"""
    diff = np.sum(l_new) - np.sum(l_ref)
    if diff < 0: #if the total degree is inferior to the one of the reference, we add uniformly at random the difference in the sequence
        l_index_sup = np.random.randint(low=0,high=len(l_new),size=abs(diff))
        for index in l_index_sup:
            l_new[index] += 1
    elif diff > 0: #if the total degree is superior to the one of the reference, we remove uniformly at random the difference in the sequence
        l_index_sup = np.random.randint(low=0, high=len(l_new), size=diff)
        for index in l_index_sup:
            l_new[index] -= 1
    else: #if already ok do nothing
        return l_new
    return l_new

def random_multidegree_sequence(path):
    """Generate nodewise/layerwise random multidegree sequence to mimic the generation of random multilayer networks
    fitting with the input multilayer network in terms of number of nodes, layers and density. Then the function computes
    the theoretical multidegree distance between the multidegree sequence and its averaged over the outputs of the stochastic rewiring algorithm
            Input:
            path -- path where the multiplex network is
            Output:
            d_x --  theoretical  nodewise distance  between the input multiplex network's multidegree sequence and its theoretical averaged counterpart
            d_y -- theoretical  layerwise distance  between the input multiplex network's multidegree sequence and its theoretical averaged counterpart
            l_dx_random -- list of theoretical layerwise distance between the input random multidegree sequences and their theoretical averaged counterpart
            l_dy_random -- list of theoretical nodewise distance between the input random multidegree sequences and their theoretical averaged counterpart
            N -- number of nodes
            M -- number of layers
            density -- desnity of the input multiplex network
            lK_x -- nodewise multidegree sequence of the input multiplex network
            lK_y -- layerwise multidegree sequence of the input multiplex network
     """
    with open(path,'rb') as f:
        parameter = pickle.load(f)
    lK_x = parameter['l_K']
    lK_y = parameter['l_Ktranspose']
    # to get an uniform probability with mplex
    l_d_x_random = []
    l_d_y_random = []
    N = parameter['N']
    M = parameter['M']
    density = sum(lK_x) / (M * N * (N - 1))
    p_node = (N * (N - 1) - 2) / (M * N * (N - 1) - 2) # p_node parameters to get an uniform probability for a link to be rewired in mplex
    p_layer = 2 * (M - 1) / ((M * N * (N - 1) - 2)) # p_layer parameters to get an uniform probability for a link to be rewired in mplex
    #print(p_layer)
    p_tel = 1 - p_node - p_layer
    #print(density)
    d_x,d_y = compute_theoretical_distance(lK_x, lK_y, 1, [p_node, p_layer, p_tel]) #theoretical distance between network and averaged version of it after rewiring
    #print(d_x,d_y)
    for u in range(number_iteration): #build random nodewise and layerwise multidegree sequence
        #print(u)
        l_K_random = []
        l_K_DS_random = []
        for i in range(M):
            k_DS_random = 2*np.random.binomial(N*(N-1)/2,density) #layerwise random multidegree of layer i
            l_K_DS_random.append(k_DS_random)
        for i in range(N):
            k_random = np.random.binomial(M * (N - 1), density) #nodewise random multidegree of layer i
            l_K_random.append(k_random)
        l_K_random = equalizer(l_K_random,lK_x) #functions to get the same exact density as the empirical network
        l_K_DS_random = equalizer(l_K_DS_random, lK_y)
        d_x_random,d_y_random = compute_theoretical_distance(l_K_random,l_K_DS_random,1,[p_node,p_layer,p_tel])
        l_d_x_random.append(d_x_random) #disatnce between generated random multidegree sequence and its averaged version after rewiring.
        l_d_y_random.append(d_y_random)
    return d_x,d_y,l_d_x_random,l_d_y_random,N,M,density,lK_x,lK_y

if __name__ == "__main__":
    ##Parameters
    path = "mplex_sparse_random_multilayer_N_200_M_200_nnz_39218.dat"  #path where the initial multiplex network is
    number_iteration = 100 #number of distances over which the average distance is made
    writing_path = "normalized_distance/csv_normalized_distance.csv" #path where results are saved.
    ##RUN
    d_x,d_y,l_dx_random,l_dy_random,N,M,density,lK_x,lK_y = random_multidegree_sequence(path)
    print(np.std(lK_y))
    ##Save data
   # Notably we save the ratio d_x/avg(l_dx_random) and d_y/avg(l_dy_random)
    with open(writing_path,'a') as f:
        string_csv = path + "," + str(N) + "," + str(M) +"," +str(d_x/np.mean(l_dx_random)) + ","+ str(d_y/np.mean(l_dy_random))+","+str(np.mean(lK_x))+","+str(np.mean(lK_y))+","+ str(density)+"\n"
        f.write(string_csv)

