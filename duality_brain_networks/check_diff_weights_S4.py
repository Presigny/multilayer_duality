"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import os
import matplotlib.pyplot as plt

def to_supra_index(N,uplet):
    """Convert tensorial indices (layer,layer,node,node) into flattened matrix indices (row,col)
            Input:
            N -- number of nodes
            uplet -- matrix indices (row,col)
            Output:
            row, col -- matrix indices (layer,layer,node,node) """
    k, l, i, j = uplet[0],uplet[1],uplet[2],uplet[3]
    row = k*N +i
    col =  l*N +j
    return row,col

##Parameters
save_path = "csv_matrix/"
number = 70 #number of layers of the considered reduced matrix
load_path = "preprocessing/database_reduced/%s_sym_individual_matrix_PAT/"%(number)
l_load_path = os.listdir(load_path) #list the patient reduced matrices
load_HC_average = "preprocessing/database_reduced/%s_sym_average_SUJ_matrix/%s_sym_average_SUJ.dat"%(number,number)
N=70 #number of nodes of the non reduced matrices
M=77 #number of layers of the non reduced matrices
##
l_intra = [] #select intralayer number of links
l_inter = [] #select interlayer number of links
l_replica = [] #select replica number of links
A_average = np.load(load_HC_average) # matrix of the average healthy subject
for matrix in l_load_path:
    A_matrix = np.load(load_path+matrix)
    B = A_matrix - A_average
    B = np.abs(B) #this matrix summarises the absolute weight difference between aptient matrix and average healthy subject matrix
    sum_intra = 0
    sum_inter = 0
    sum_replica = 0
    for k in range(M):
        for l in range(M):
            row_entrance, col_entrance = to_supra_index(N,[k,l,0,0]) #first coordinates in an (inter)layer
            row_sortance,col_sortance = to_supra_index(N,[k,l,N,N]) #last coordinates in an (inter)layer
            if k==l:
                sum_intra += np.sum(B[row_entrance:row_sortance,col_entrance:col_sortance]) #sum the weights contained in a layer
            else:
                tr = np.trace(B[row_entrance:row_sortance, col_entrance:col_sortance])
                sum_inter += np.sum(B[row_entrance:row_sortance,col_entrance:col_sortance])-tr #sum the weights that containe in an interlayer, subtratcing the one from replica links
                sum_replica += tr #sum the weights that links replica nodes
    l_intra.append(sum_intra)
    l_inter.append(sum_inter)
    l_replica.append(sum_replica)
    print(sum_intra,sum_inter,sum_replica)
    print(np.sum(B)) #difference in the total wieghts between the considered matrix and the matrix of the average healthy subject
    #np.savetxt(save_path+matrix+".csv",A_matrix,delimiter=',')
    print(matrix+" DONE")

##Make the histogramm of intralayer, interlayer nad replica weights
alpha= 0.5
plt.hist(l_intra,label="intralayer",color="blue",alpha=alpha)
plt.hist(l_inter,label="interlayer",color="orange",alpha=alpha)
plt.hist(l_replica,label="replica",color="green",alpha=alpha)
plt.xscale("log")
#plt.title(str(number)+" layers")
plt.legend()
plt.show()
