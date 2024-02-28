"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import pickle
import numpy as np
###################################################################
# Description: Take the node and layer multilayer degree of all multilayer
#               networks in a given directory and compute and store
#               their standard deviations and the density of the network


###################################################################
#l_path = ["mplex_sparse_random_multilayer_N_200_M_200_nnz_39218.dat"]
l_path = ["preprocessing/database/Arxiv/arxiv_multiplex.dat"]
writing_path = "standard_deviation.csv" #select the file in which you write
type = "multiplex" #select the type of multilayer

with open(writing_path, 'w') as f: #intialize the file with title of each column
    string_csv = "name,N,M,std_x,std_y,density\n"
    f.write(string_csv)

##RUN

for path in l_path:
    with open(path,'rb') as f: #Load the data
        parameter = pickle.load(f)
    A = parameter['sparse']
    lK_x = parameter['l_K']
    lK_y = parameter['l_Ktranspose']
    N = parameter['N']
    M = parameter['M']
    #Compute the standard deviations
    std_x = np.std(lK_x)
    std_y = np.std(lK_y)
    if type == "multiplex":
        density = 2 * A.nnz / (M * N * (N + M - 2))
        #Save the data
        # with open(writing_path,'a') as f: #save name,N,M,std_X,std_Y,density
        #     string_csv = path + "," + str(N) + "," + str(M) +"," +str(std_x) + ","+ str(std_y)+","+ str(density)+"\n"
        #     f.write(string_csv)



