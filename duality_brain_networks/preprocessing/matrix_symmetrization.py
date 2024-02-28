"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import os

def symmetrize_matrix(A):
    """
    Symmetrize an adjacency matrix by replacing the value of the weakest link to the strongest one
    Input:
    A -- adjacency matrix
    Output:
    A -- symmetrized adjacency matrix
    """
    for i in range(len(A)):
        for j in range(len(A[0])):
            if A[i,j] > A[j,i]:
                A[j,i] = A[i,j]
            else:
                A[i, j] = A[j,i]
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

# load_path = "database/individual_matrix/"
# save_path = "database/sym_individual_matrix/"

# l_load_path = os.listdir(load_path)
# for group in l_load_path:
#     A = np.load(load_path +group)
#     A = symmetrize_matrix(A)
#     #Save the data
#     with open(save_path +group+"_sym", 'wb') as f:
#         np.save(f,A)
#     print(group," DONE")
## USe this for Fig S4
number = 70
load_path = "database_reduced/%s_sym_individual_matrix/"%(number)
save_path = "database_reduced/%s_sym_individual_matrix_PAT/"%(number)
l_load_suj,l_load_pat = group_by_type(load_path)
for matrix in l_load_pat:
    A = np.load(load_path + matrix)
    A = symmetrize_matrix(A)
    # Save the data
    with open(save_path + matrix + "_sym", 'wb') as f:
        np.save(f, A)
    print(matrix, " DONE")
