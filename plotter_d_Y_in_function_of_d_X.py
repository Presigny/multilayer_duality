"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib.collections import LineCollection
##################################################
# in readme
# We computed 100 couple of multidegree distances (d_X and d_Y) from 100 synthetic random multilayer networks.
# Each couple of multidegree distance is associated with a different R (proportion of links rewired).
# Therefore, each synthetic random multilayer network is associated with a different R.
# The previous process was applied for 3 different combinations of p_node, p_layer, p_tel
##################################################
##Path for data
#path = 'figure_2_full_model_2_N_200_M_200_stat.dat'
path = 'mock_data/figure_2_full_model_2_N_200_M_200_25.dat'#figure_2_full_model_2_N_200_M_200.dat' #full multilayer case
#path = 'figure_2_mplex_model_2_N_200_M_200_stat.dat' #multiplex case
with open(path, 'rb') as f:
    parameter = pickle.load(f) #load the data's dictionary
## Load data
l_dD0 = parameter['l_distance_K0'] # list of theoretical distances (d_X) associated with each synthetic random multilayer network for the 3 combinations [p_node,p_layer,p_tel]
l_dDS = parameter['l_distance_KDS'] # list of theoretical distances (d_Y) associated with each synthetic random multilayer network for the 3 combinations [p_node,p_layer,p_tel]
l_proba_event = parameter['l_proba_event'] # combinations  [p_node,p_layer,p_tel]
N = parameter['N'] #number of nodes
M = parameter['M'] #number of layers
l_R = parameter['l_R'] #list of proportion of links that were rewired
l_density = parameter['density'] #list of densities of the networks
l_dD0_th = np.zeros([len(l_proba_event),len(l_R)]) # list of average theoretical distance for random multilayer networks (d_X) for 3 combinations [p_node,p_layer,p_tel]
l_dDS_th = np.zeros([len(l_proba_event),len(l_R)]) # list of average theoretical distance for random multilayer networks (d_Y) for 3 combinations [p_node,p_layer,p_tel]
l_color = ["red","blue","orange"]
# # FOR FULL MULTILAYER
# computing with formula the average theoretical distances for random multilayer networks for 3 combinations [p_node,p_layer,p_tel]
for i in range(len(l_proba_event)):
    for k in range(len(l_R)):
        l_dD0_th[i][k] = N*M*l_R[k]*(1-l_proba_event[i][0])*np.sqrt(l_density[i][k]*(1-l_density[i][k]))
        l_dDS_th[i][k] = M *N* l_R[k] * (1 - l_proba_event[i][1])*np.sqrt(l_density[i][k]*(1-l_density[i][k]))

## FOR MULTIPLEX NETWORKS, uncomment what is below and comment what is above
# computing with formula the average theoretical distances for random multiplex networks for 3 combinations [p_node,p_layer,p_tel]
# for i in range(len(l_proba_event)):
#     for k in range(len(l_R)):
#         l_dD0_th[i][k] = N*np.sqrt(M)*l_R[k]*(1-l_proba_event[i][0])*np.sqrt(l_density[k]*(1-l_density[k]))
#         l_dDS_th[i][k] = np.sqrt(2)*N * np.sqrt(M) * l_R[k] * (1 - l_proba_event[i][1])*np.sqrt(l_density[k]*(1-l_density[k]))

##PLOT
fig, ax = plt.subplots(figsize=(6,6.9))
plt.rcParams['font.size'] = '18'
for i in range(len(l_proba_event)): #loop over the 3 combinations [p_node,p_layer,p_tel]
    plt.scatter(l_dD0[i],l_dDS[i],s=10,c=l_R,cmap='viridis')
    points = np.array([l_dD0_th[i], l_dDS_th[i]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(np.min(l_R), np.max(l_R))
    lc = LineCollection(segments, cmap='viridis',norm=norm) #build the average theoretical lines
    lc.set_array(l_R)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
fig.colorbar(line)
fontsize = '18'
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
ax.set_xlim(np.min(l_dD0_th), np.max(l_dD0)+5) #+5 to free some space above
ax.set_ylim(np.min(l_dDS_th), np.max(l_dD0)+5)
#ax.set_aspect(1. / ax.get_data_ratio(), adjustable='box')
ax.spines[['right', 'top']].set_visible(False)
#plt.savefig("2b_mplex.png",dpi=300,format='png') #save the figure
plt.show()
