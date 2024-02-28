import numpy as np
import os
import pickle
import pandas as pd
from sklearn.preprocessing import normalize

def check_non_normalized_entry(A):
    for i in range(len(A)):
        if max(A[i]) > 1:
            return True

save_path = "csv_matrix/"
load_path = "Database/individual_matrix/"
l_load_path = os.listdir(load_path)

for matrix in l_load_path:
    A_matrix = np.load(load_path+matrix)
    np.savetxt(save_path+matrix+".csv",A_matrix,delimiter=',')
    print(matrix+" DONE")

