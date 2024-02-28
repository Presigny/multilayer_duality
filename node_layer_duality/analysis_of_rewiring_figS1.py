"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import pickle
import os
import numpy as np
from math import fsum
import matplotlib.pyplot as plt

def from_file_curation(path):
    """Data curation to supress the skipline, the [, the ] in the files
        Input:
        path -- path of the file to be read
        Output:
        m_K -- Array of data extracted from the files"""
    m_K = []
    file1 = open(path, 'r')
    Lines = file1.readlines()
    file1.close()
    for i in range(len(Lines)):
        #replace useless characters
        Lines[i] = Lines[i].replace('\n', '')
        Lines[i] = Lines[i].replace('[', '')
        Lines[i] = Lines[i].replace(']', '')
        Lines[i] = Lines[i].split(sep=",")
        m_K.append([int(Lines[i][j]) for j in range(len(Lines[i]))])
    return m_K

def get_meta_data(path):
    """Read the name of file to obtain metadata on it
        Input:
        path -- path of the file to be read
        Output:
        rule -- index of the rule that has been used
        N -- number of nodes
        M -- number of layers
        R -- proportion of links that has been rewired
        density -- density of the multilayer networks"""
    split_path = path.split(sep="_")
    rule = split_path[4]
    N = int(split_path[6])
    M = int(split_path[8])
    R = float(split_path[10])
    density = float(split_path[12])
    return rule,N,M,R,density

def distance(l_K1,l_K2):
    """Return the Euclidean distance between l_K1 and l_K2
            Input:
            l_K1 -- multidegree vector
            l_K2 -- multidegree vector
            Output:
            Euclidean distance between l_K1 and l_K2"""
    l_K  = np.subtract(l_K1,l_K2)
    return np.sqrt(np.dot(l_K,l_K))

def load_initial_multilayer(which,N,M,number_of_edge):
    """
    Return the initial multidegree sequences from which all the rewiring sequences in the script were done
    :param which: if "Nodewise", select the node multidegree sequence, if "Layerwise, select the layer multidegree sequence
    :param N: number of nodes
    :param M: number of layers
    :param number_of_edge: number of links
    :return The layer or node multidegree sequence on the initial multilyaer network
    """
    name_path = "full_lil_sparse_random_multilayer_N_%s_M_%s_nnz_%s_inter_True_mplex_random_True_no_iso_2.dat" % (
    N, M, number_of_edge) #path have to follow this nomenclature for the script to work
    with open(name_path, 'rb') as f:
        parameter = pickle.load(f)  # we load the data file containing almost all the parameters already
    # parameters characterizing the random multilayer networks
    multi_degree = parameter['l_K']
    multi_degree_DS = parameter['l_Ktranspose']
    if which == "Layerwise":
        return multi_degree_DS
    if which == "Nodewise":
        return multi_degree

def theoretical_result(which,initial,idx,M,N,number_edges,n):
    """

    :param which: if "Nodewise", select the node multidegree sequence, if "Layerwise, select the layer multidegree sequence
    :param initial: initial multidegree sequence (depends on "which")
    :param idx: fix the rewiring parameters
    :param M: number of layers
    :param N: number of nodes
    :param number_edges: number of links
    :param n: number of links to rewire
    :return: theoretical distance for full multilayer network
    """
    number_pair = M * M * N * (N - 1) / 2  # full multilayer case
    pair_links = M * M  # full multilayer case
    pair_layers = (N * (N - 1)) / 2
    pair_tel = number_pair - pair_layers - pair_links + 1
    if idx == 0:
        p_space=1
        p_edge=0
        p_tel=0

    elif idx==64:
        p_space = 0
        p_edge = 1
        p_tel = 0

    elif idx==74:
        p_space = 0
        p_edge = 0
        p_tel = 1

    if which == "Layerwise":
        l_K_C = [0] * M
        for k in range(M):
            l_K_C[k] = initial[k] + n * p_edge * (
                    2 * M / (pair_links - 1) - ((pair_links) * initial[k]) / (number_edges * (pair_links - 1)))
            l_K_C[k] += (n / (1 - number_edges / number_pair)) * p_tel * (
                2 * M * (pair_layers - 1) / pair_tel - initial[k] * (number_pair - pair_links) / (
                pair_tel * number_edges))
    else:
        l_K_C = [0] * N
        for i in range(N):
            l_K_C[i] = initial[i] + n * p_space * (
                    (N - 1) / (pair_layers - 1) - ((pair_layers) * initial[i]) / (
                    number_edges * (pair_layers - 1)))
            l_K_C[i] += n * p_tel * (
                    (N - 1) * (pair_links - 1) / pair_tel - initial[i] * (number_pair - pair_layers) / (
                    pair_tel * number_edges))

    return distance(l_K_C,initial)

def rule_selector(idx):
    if idx==0:
        return r"$p_{space}=1$"
    elif idx==64:
        return r"$p_{edge}=1$"
    elif idx==74:
        return r"$p_{tel}=1$"


####################################################################################################
"""Parameters to display the  theoretical and numerical distance in function of the rewiring parameter (idx) and 
 the proportion of links that are rewired (l_R), expected sequences are computed with a least 100 realizations"""
which_description = "Layerwise"#The description to display in the final plot (Nodewise or Layerwise)
l_R = []
l_density = []
N = 20
M = 20
l_K_avg = [0]*M #change the cardinality here: if node, put N; if layer put M
l_D = []
l_D_uncertainty = []
l_D_th = []
idx = 64 #select the type of rewiring
path = "mock_data/sparse_rewiring_data/" #path where the data are


####################################################################################################
l_restore_path = os.listdir(path)
l_restore_path = [l_restore_path[i] for i in range(len(l_restore_path)) if l_restore_path[i].count(which_description)==1
                  and l_restore_path[i].count("_N_%s"%(str(N)))==1 and l_restore_path[i].count("_M_%s"%(str(M)))==1
                  and l_restore_path[i].count("_rule_%s"%(str(idx)))==1 ] #select precisely the data that fit the above parameters
for path_distribution in l_restore_path: # for each data file in the path
    rule,N,M,R,density = get_meta_data(path+path_distribution) #extract the data from the files
    number_of_edge = int(density*(N*M*(N*M-1))/2)
    l_R.append(R)
    l_density.append(density)
    m_K = np.array(from_file_curation(path+path_distribution))
    m_K_t = m_K.transpose()
    sample_size = len(m_K) #size of the number of computed multidegree sequence present in path_distribution (for a given condition)
    print(sample_size)
    for j in range(sample_size):
        l_K_avg = [l_K_avg[i]+m_K_t[i][j] for i in range(len(l_K_avg))] #expected multidegree sequence over the sample size
    l_K_avg = [l_K_avg[i]/sample_size for i in range(len(l_K_avg))]
    SEM_K = [np.std(m_K_t[i])/np.sqrt(sample_size) for i in range(len(l_K_avg))] #standard error of mean (SEM-
    Uncertainty_K = [l_K_avg[i]+3*SEM_K[i] for i in range(len(l_K_avg))] #uncertainty from the SEM
    initial_K = load_initial_multilayer(which_description, N, M, number_of_edge)
    print("sum K_avg: ", fsum(l_K_avg))
    print("sum K_initial ", fsum(initial_K))
    D = distance(l_K_avg,initial_K) # distance between expected multidegree sequence and initial one
    D_th = theoretical_result(which_description,initial_K,idx,M,N,number_of_edge,number_of_edge*R)
    l_D.append(D) # numerical distance for a given R
    l_D_th.append(D_th) #theoretical distance for a given R
    Uncertainty_D = distance(Uncertainty_K,initial_K) #uncertainty on the distance
    l_D_uncertainty.append(abs(D-Uncertainty_D))
    if which_description == "Layerwise":
        l_K_avg = [0] * M
    else:
        l_K_avg = [0] * N

####reorder####

set_density = set(l_density)
r_density = list(set_density)
set_density = np.sort(r_density)
set_R = set(l_R)
dico_R = {}
dico_distance = {}
dico_distance_th = {}
dico_uncertainty = {}

for density in set_density:
    dico_R[str(density)] = []
    dico_distance[str(density)] = []
    dico_uncertainty[str(density)] = []
    dico_distance_th[str(density)] = []

for i in range(len(l_restore_path)):
    dico_R[str(l_density[i])].append(l_R[i])
    dico_distance[str(l_density[i])].append(l_D[i])
    dico_uncertainty[str(l_density[i])].append(l_D_uncertainty[i])
    dico_distance_th[str(l_density[i])].append(l_D_th[i])
### Plot the theoretical and numerical distance in fucntion of the proportion of links to be rewire
alpha=0.5
u = 0
color= ['blue','orange','green','red','purple']
displayed_density = [1e-3, 5e-3, 0.01, 0.05, 0.1]
plt.rcParams['font.size'] = '26'
fig,ax = plt.subplots(figsize=(11,11))
#ax.spines[['right', 'top']].set_visible(False)
for density in set_density:
    plt.plot(dico_R[str(density)], dico_distance_th[str(density)],color=color[u],alpha=0.7)
    plt.scatter(dico_R[str(density)], dico_distance[str(density)], label=r"$\rho=$" + str(displayed_density[u]),color=color[u],alpha=alpha,s=100)
    plt.errorbar(dico_R[str(density)], dico_distance[str(density)], yerr=dico_uncertainty[str(density)],fmt='o',color=color[u],alpha=alpha)
    u += 1 # for the color to change
plt.xlabel(r"$r$",fontsize='32')
if which_description == "Nodewise":
    plt.ylabel(r"$d_\mathcal{X}$",rotation=90,fontsize='32',labelpad=10)
else:
    plt.ylabel(r"$d_\mathcal{Y}$", rotation=90, fontsize='32', labelpad=10)

rule = rule_selector(idx)
title = "Full multilayer Erdös-Renyi N=%s, M=%s "%(str(N),str(M))+rule_selector(idx)+" "+which_description
plt.legend()
ax = plt.gca()
ax.set_aspect(1/ax.get_data_ratio(),adjustable='box')
plt.savefig(title,transparent="True")
plt.show()

