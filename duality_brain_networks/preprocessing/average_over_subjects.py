"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import os

def averaging_array(load_path,l_group):
    """Compute the average supra-adjacency matrix of a given set of supra_adjacency_matrix
        Input:
        load_path -- path of the directory that contains the group
        l_group -- filemanes of the considered group
        Output:
        A -- average supra-adjacency matrix"""
    A = np.load(load_path + l_group[0])
    a = A[0][0]
    for i in range(1,len(l_group)):
        B = np.load(load_path + l_group[i])
        A = np.add(A,B)
        a += B[0][0]
    A = A/len(l_group)
    a = a/len(l_group)
    print(a)
    return A

def group_by_type(load_path):
    """Separe the filenames of Patients and Subjects
        Input:
        load_path -- global path where patients and subjects names are mixed
        Output:
        l_group_suj -- subjects filename
        l_group_pat -- patients filename"""
    l_load_path = os.listdir(load_path)
    l_group_suj = []
    l_group_pat = []
    for name_1 in l_load_path:
        if name_1.count("PP") != 0:#detect patients filename
            l_group_pat.append(name_1)
        else:
            l_group_suj.append(name_1)
    return l_group_suj,l_group_pat


load_path = "database_reduced/70_sym_individual_matrix/"
save_path = "database_reduced/70_sym_average_SUJ_matrix/70_sym_average_SUJ.dat"


#load_path = "database/sym_individual_matrix/"
#save_path = "database/sym_avgSUJ_matrix/sym_average_SUJ.dat"

l_load_suj,l_load_pat = group_by_type(load_path)
A = averaging_array(load_path,l_load_suj)
#Save the average supra-adjacency matrix
with open(save_path, 'wb') as f:
    np.save(f,A)
print(" DONE")