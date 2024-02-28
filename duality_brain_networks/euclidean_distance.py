"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import os
import pickle
import numpy as np

def compute_distance(l1,l2):
    delta = np.subtract(l1, l2)
    return np.sqrt(np.dot(delta, delta))

## Uncomment if you want to compute euclidean distance between patients and average subjects with th reduced matrix
#number_layer = 7
#main_path = "Database_reduced/"
#patient_path = main_path+"%s_sym_multistrength/"%(number_layer)
#load_path_avg = main_path+"multistrength_%s_sym_average_HC.dat"%(number_layer)
#writing_path = main_path+"%s_std_to_avg_HC_sym.csv"%(number_layer)

## euclidean distance between patients and average subjects
main_path = "preprocessing/database_reduced/"
patient_path = main_path+"70_sym_multistrength_PAT/" #select the file where the multistrength of patients are
load_path_avg = main_path+"/70_sym_multistrength_avgSUJ/70_sym_average_SUJ.dat_multistrength" #select the file where the multistrength of the average healthy subject are
writing_path = main_path+"70_euclidean_distance_to_avg_SUJ_sym.csv" #where the data will be saved

l_restore_path = os.listdir(patient_path)
l_name = [l_restore_path[i][:-18] for i in range(len(l_restore_path))]
data = {}
data_DS = {}
with open(writing_path, 'w') as f:  # save name,N,M,std_X,std_Y,density
        string_csv = "name,d_X,d_Y\n"
        f.write(string_csv)

with open(load_path_avg, 'rb') as f:  # load average multistrength
    dico = pickle.load(f)
multistrength_avgHC = dico['multistrength']
multistrength_DS_avgHC = dico['multistrength_DS']

for strength_path_2 in l_restore_path: #for every patients
    container = []
    container_DS = []
    with open(patient_path + strength_path_2, 'rb') as f: # load average multistrength
        dico = pickle.load(f)
    multistrength = dico['multistrength']
    multistrength_DS = dico['multistrength_DS']

    distance = compute_distance(multistrength_avgHC,multistrength) #compute euclidean distance between average healthy and patients
    distance_DS = compute_distance(multistrength_DS_avgHC,multistrength_DS)
    with open(writing_path, 'a') as f:  # save name,N,M,std_X,std_Y,density
        string_csv = str(strength_path_2[:-18]) + "," + str(distance) + "," + str(distance_DS) + "\n"
        f.write(string_csv)
