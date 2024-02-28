"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import ranksums

if __name__ == "__main__":
    ## Parameters
    #l_element= [3,7,15,19,23,31,38,46,54,61,70] #list of the number of layers
    l_element = [70]  # list of the number of layers
    l_load_path = []
    l_mean_X = [] #mean of nodewise distance between subjects and patients
    l_SEM_X = [] #standard error of mean (nodewise)
    l_mean_Y = [] #mean of layerwise distance between subjects and patients
    l_SEM_Y = [] #standard error of mean (layerwise)
    l_p_value = [] # ranksum between nodewise distance distribution and layerwise one
    ## Load the euclidean distance obtained for each group of reduced matrices
    for number in l_element:
        load_path = "preprocessing/database_reduced/%s_euclidean_distance_to_avg_SUJ_sym.csv" %number
        df = pd.read_csv(load_path)
        distance_suj_pat = df['d_X'].values.tolist()
        distance_suj_pat_DS = df['d_Y'].values.tolist()
        l_mean_X.append(np.mean(distance_suj_pat))
        l_mean_Y.append(np.mean(distance_suj_pat_DS))
        l_SEM_Y.append(np.std(distance_suj_pat_DS) / np.sqrt(len(distance_suj_pat))) #standard error of mean
        l_SEM_X.append(np.std(distance_suj_pat) / np.sqrt(len(distance_suj_pat)))
        print(number,"ranksum ",ranksums(distance_suj_pat, distance_suj_pat_DS)[1]) #use the ranksums test for significance between node and layer dimension
        l_p_value.append(ranksums(distance_suj_pat, distance_suj_pat_DS)[1])
    ## Load the original 77 layers networks
    l_element.append(77)
    load_path = "preprocessing/database/euclidean_distance_to_avg_HC_sym.csv"
    df = pd.read_csv(load_path)
    distance_suj_pat = df['d_X'].values.tolist()
    distance_suj_pat_DS = df['d_Y'].values.tolist()
    l_mean_X.append(np.mean(distance_suj_pat))
    l_mean_Y.append(np.mean(distance_suj_pat_DS))
    l_SEM_Y.append(np.std(distance_suj_pat_DS) / np.sqrt(len(distance_suj_pat)))
    l_SEM_X.append(np.std(distance_suj_pat) / np.sqrt(len(distance_suj_pat)))
    print(77,"ranksum ",ranksums(distance_suj_pat, distance_suj_pat_DS)[1]) #use the ranksums test for significance between node and layer dimension
    l_p_value.append(ranksums(distance_suj_pat, distance_suj_pat_DS)[1])

    ## make the plot mean distances in function of M
    fontsize= '20'
    fig,ax = plt.subplots(figsize=(5,5))
    ax.spines[['right', 'top']].set_visible(False)
    plt.yscale("log")
    plt.rcParams.update({"font.size":20})
    plt.errorbar(l_element, l_mean_X, yerr=l_SEM_X, xerr=None, color ='#0047AB')
    plt.plot(l_element, l_mean_X, color='#BDD7EE')
    plt.scatter(l_element, l_mean_X, color='#BDD7EE', label=r"$\langle d_X \rangle$", s=100, edgecolors='black')
    plt.errorbar(l_element, l_mean_Y, yerr=l_SEM_Y, xerr=None, color ='#CC5500')
    plt.plot(l_element, l_mean_Y, color='#F8CBAD')
    plt.scatter(l_element, l_mean_Y, color='#F8CBAD', label=r"$\langle d_Y \rangle$", s=100, edgecolors='black')
    plt.xlabel("M",fontsize=fontsize)
    plt.ylabel("Distance",rotation=90,fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.tight_layout()
    plt.legend()
    ax = plt.gca()
    ax.set_aspect(1/ax.get_data_ratio(),adjustable='box')
    plt.savefig("distance_func_M_euclidean",transparent=True)
    plt.show()
