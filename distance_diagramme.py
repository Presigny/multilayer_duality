"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""

import matplotlib.pyplot as plt
import pickle
from numpy import mgrid
from math import fsum
from scipy.interpolate import griddata
import numpy as np
import matplotlib as mpl
if __name__ == "__main__":
    ##Parameter of the plot
    data_path = 'mock_data/distance_diagramme/a_full_model_2_N_200_M_200_R_39619.0_edges_79238_pair_links_40000.dat'
    ##Load data
    with open(data_path, 'rb') as f:
        data = pickle.load(f)
    N = data['N']
    M = data['M']
    R = data['R']
    l_distance_K0 = data['l_distance_K0']
    l_distance_KDS = data['l_distance_KDS']
    l_distance_K0 = np.array(l_distance_K0) #nodewise distance for each tested probability combinations
    l_distance_KDS = np.array(l_distance_KDS) #layerwise distance for each tested probability combinations
    l_proba_event = data['l_proba_event']
    l_proba_event = np.array(l_proba_event)

    p_intra = [uplet[0] for uplet in l_proba_event] # extract p_node parameters
    p_inter = [uplet[1] for uplet in l_proba_event] # extract p_layer parameters
    p_tel = [uplet[2] for uplet in l_proba_event] # extract p_tel parameters

    points = [[p_intra[i],p_inter[i]] for i in range(len(p_intra))]
    ##Build the grid
    grid_x, grid_y = mgrid[0:1:100j, 0:1:100j]
    grid_z1 = griddata(points, values=l_distance_K0, xi=(grid_x, grid_y), method='linear')
    grid_z2 = griddata(points, values=l_distance_KDS, xi=(grid_x, grid_y), method='linear')

    if max(l_distance_KDS)-max(l_distance_K0)< 0:
        upper_bound = max(l_distance_K0)
    else:
        upper_bound = max(l_distance_KDS)
    my_cmap = plt.get_cmap('Spectral_r') #set the colormap
    my_cmap.set_bad('w')
    fig, ax = plt.subplots(figsize=(20,6),sharey='row')
    plt.rcParams['font.size'] = '24'
    ## Subplot of the nodewise distance diagramme
    plt.subplot(121)
    my_cmap = plt.get_cmap('Spectral')
    plt.imshow(grid_z1.T, extent=(0, 1, 0, 1), origin='lower', cmap=my_cmap)
    plt.title("$d_\mathcal{X}$",y=1.1)
    plt.xlabel("$p_{node}$",fontsize=24)
    plt.ylabel("$p_{layer}$", rotation=90, labelpad=20, fontsize=24)
    for i in range(1,11):
       plt.axline([0,i/10],[i/10,0], color='black',linewidth=0.4)
    plt.colorbar(pad=0.1)
    plt.clim(-upper_bound,upper_bound)
    plt.xticks([0.2, 0.4, 0.6, 0.8, 1])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    #plt.show()
    ## Subplot of the layerwise distance diagramme
    plt.subplot(122)
    my_cmap = plt.get_cmap('Spectral_r')
    plt.imshow(grid_z2.T, extent=(0, 1, 0, 1), origin='lower', cmap=my_cmap)
    plt.title("$d_\mathcal{Y}$",y=1.1)
    plt.xlabel("$p_{node}$",fontsize=24)
    for i in range(1,11):
       plt.axline([0,i/10],[i/10,0], color='black',linewidth=0.4)
    plt.colorbar(pad=0.1)
    plt.xticks([0.2, 0.4, 0.6, 0.8, 1])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.clim(-upper_bound,upper_bound)
    plt.tight_layout(pad=1.5)
    #plt.savefig("distance_diagram",transparent=True)
    plt.show()

