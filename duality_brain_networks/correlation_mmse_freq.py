"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests
import mne
from scipy.sparse import csr_matrix


def clean_scores(dataframe):
    """Transform all the ND indices in the dataframes into -1
    Input:
    dataframe -- pandas dataframe
    Output:
    dataframe -- cleaned pandas dataframe"""
    for ind in dataframe.index:
        for col in dataframe.columns:
            if dataframe.loc[ind, col] == 'nd' or dataframe.loc[ind, col] == 'ND':
                dataframe.loc[ind, col] = -1
    return dataframe

def score_in_df(dataframe_1, dataframe_2, score):
    """For a given score, insert its related value for each patient from dataframe_2 to dataframe_1
    Input:
    dataframe_1 -- pandas dataframe
    dataframe_2 -- pandas dataframe
    score -- congitive score considered
    Output:
    dataframe_1 -- pandas dataframe with score associated to each patient"""
    dataframe_1.insert(len(dataframe_1.columns), str(score), " ")
    for ind in dataframe_1.index:
        name = re.findall(r'\d+', ind)
        if len(name[0]) > 2:
            name[0] = name[0][1:]
        a = dataframe_2[dataframe_2['Code_IRM'].str.contains(name[0])]
        if len(a) != 0:
            dataframe_1.at[ind, score] = int(a[score])
        else:
            dataframe_1.at[ind, score] = -1
    return dataframe_1

def recalibrate(df, name_score):
    """Transform the dataframe into 3  associated lists
        Input:
        df -- pandas dataframe
        name_score -- name of the considered cognitive score
        Output:
        new_dist_X -- nodewise multidegree of ROI, ROI_index
        new_dist_Y -- layerwise multidegree of frequency, layer_index
        new_score -- cognitive score of patients"""
    l_dist_X = df['K_X_' + str(ROI_index)].values.tolist()
    l_dist_Y = df['K_Y_' + str(layer_index)].values.tolist()
    l_score = df[name_score].values.tolist()
    new_dist_X = []
    new_dist_Y = []
    new_score = []
    for dist_X, dist_Y, score in zip(l_dist_X, l_dist_Y, l_score):
        if score != -1:
            new_dist_X.append(dist_X)
            new_dist_Y.append(dist_Y)
            new_score.append(score)
    return new_dist_X, new_dist_Y, new_score


def remove_outliers(l_X, l_Y, l_sco):
    """Remove the outliers with abs(z-score) < 3 and update the node/layerwise multidegree distributions
    and the score list and provide the index and the zscore of the outliers
    """
    z_dx = zscore(l_X)
    z_dy = zscore(l_Y)
    l_outlier_x = []
    l_outlier_y = []
    for i in range(len(z_dx)):
        if abs(z_dx[i]) > 3:
            l_outlier_x.append((i, z_dx[i]))
    for i in range(len(z_dy)):
        if abs(z_dy[i]) > 3:
            l_outlier_y.append((i, z_dy[i]))
    for i in range(len(l_outlier_x)):
        l_X.pop(l_outlier_x[i][0])
        l_sco.pop(l_outlier_x[i][0])
    for i in range(len(l_outlier_y)):
        l_Y.pop(l_outlier_y[i][0])
    print(l_outlier_x, l_outlier_y)
    return l_X, l_Y, l_sco, l_outlier_x, l_outlier_y

def perm_t_stat(*args):
    """Allow to choose the Spearman correlation as a statistical function for the cluster-based permutation test
        Output:
        l_tvalue -- lsit of Spearman coefficients"""
    m_pat = args[0]
    m_score = args[1]
    l_score = m_score[:,0]
    #print(l_score)
    l_tvalue = np.zeros([len(m_pat[0])])
    l_pvalue = np.zeros([len(m_pat[0])])
    for i in range(len(m_pat[0])):
        t_values, pv = spearmanr(m_pat[:,i],l_score)
        l_tvalue[i] = t_values
        l_pvalue[i] = pv
    #t = mne.stats.f_oneway(*args)
    #print(len(t))
    return l_tvalue

def distance_freq(pos):
    """Distance between each frequency (in Hz)
    Input:
    pos -- list of frequency of each layer
    Output -- matrix of the relative position of each frequency """
    A = np.zeros([len(pos),len(pos)])
    for i in range(len(pos)):
        for j in range(len(pos)):
            A[i,j] = np.abs(pos[i]-pos[j])
    return A

l_p_value = [] #store the p values
l_spearman_coeff = [] #store the spearman coefficient
l_layer_index = range(77) # layer index of each 77 layers
matrix_pat = []
matrix_score = []
for layer_index in l_layer_index: # iterate on every layer
    ROI_index = 56 #select a random ROI_index without utility in this algorithm

    writing_path = "preprocessing/database/layer_activity_PAT.csv"
    restore_path = "preprocessing/database/sym_multistrength_PAT/"
    save_path = "preprocessing/database/plot_freq_correlation_MMSE/"
    l_restore_path = os.listdir(restore_path)

    l_name = [l_restore_path[i][26:-18] for i in range(len(l_restore_path))] #give a consistant name to generated files
    with open(writing_path, 'w') as f:  # save name,N,M,std_X,std_Y,density
        string_csv = "name,K_X_" + str(ROI_index) + ",K_Y_" + str(layer_index) + "\n"
        f.write(string_csv)
    for strength_path_2 in l_restore_path:
        with open(restore_path + strength_path_2, 'rb') as f:  # open file 1
            dic = pickle.load(f)
        multistrength = dic['multistrength'][ROI_index]
        multistrength_DS = dic['multistrength_DS'][layer_index]
        with open(writing_path, 'a') as f:  # save name,N,M,std_X,std_Y,density
            string_csv = str(strength_path_2[26:-18]) + "," + str(multistrength) + "," + str(multistrength_DS) + "\n"
            f.write(string_csv)

    #Load the excel file containing the MMSE score of patients
    df_score = pd.read_excel("metadata/PAT_metadata.xlsx")
    print(df_score.columns)
    df_activity = pd.read_csv(writing_path, encoding='unicode_escape', index_col='name')
    # df_score = df_score[['groupe','MMSE','RI\n/48','RL\n/48','R.\nIMM.\n/16','RL DIFF\n/16','RI DIFF\n/16','RT\n(L+I)','Code_IRM','âge','Sexe']]
    df_score = df_score[['groupe', 'MMSE','Code_IRM', 'âge', 'Sexe']]
    for i in df_score.index: #change the name
        df_score['Code_IRM'] = df_score['Code_IRM'].replace(df_score['Code_IRM'][i], df_score['Code_IRM'][i][17:])
    #Separate patients from subjects
    df_score_pat = df_score[df_score['Code_IRM'].str.contains('AT')]
    df_score_suj = df_score[df_score['Code_IRM'].str.contains('SUJ|Sujet')]
    df_score_pat = clean_scores(df_score_pat)
    df_activity = score_in_df(df_activity,df_score_pat,'MMSE')
    score = 'MMSE'
    l_KX,l_KY,l_score = recalibrate(df_activity,score)
    matrix_pat.append(l_KX)
    matrix_score.append(l_score)
    #l_dX,l_dY,l_score,l_outlier_X,l_outlier_Y = remove_outliers(l_dX,l_dY,l_score)
##PLOT
    fontsize=22
    plt.figure(figsize=(5,5))
    correlation = spearmanr(l_KY,l_score)
    coeff = round(correlation[0], 2)
    pvalue = round(correlation[1],4)
    l_frequency = [2 +i*0.5 for i in range(77)]
    label=r"$S$ = %s , $p$ = %s, layer index = %s Hz"%(coeff,pvalue,l_frequency[layer_index])
    plt.scatter(l_KY,l_score,color='#F8CBAD',edgecolor='black',s=100,alpha=1,label=label)
    plt.xlabel(r"$K_Y^{AD}$",fontsize=fontsize)
    plt.ylabel("MMSE",rotation=90,fontsize=fontsize)
    #plt.title(str(score)+" in function of layer multidegree disruption index")
    plt.legend()
    plt.grid()
    ax = plt.gca()
    ax.set_aspect(1/ax.get_data_ratio(),adjustable='box')
    plt.savefig(save_path+"frequency_"+str(l_frequency[layer_index])+".png",format='png')
    plt.close()
    # save the correlation between MMSE score and the total activity
    l_spearman_coeff.append(correlation[0])
    l_p_value.append(correlation[1])

    print(l_frequency[layer_index],correlation)


## Do the FDR correction on the list of p value for each frequency to correct for multiple comparisons
a,l_p_value_1,c,d = multipletests(l_p_value,alpha=0.05,method='fdr_bh')
writing_path_FDR = "preprocessing/database/FDR_0.05_freq_correlation_MMSE.csv"
with open(writing_path_FDR,"w") as f:
    for i in range (len(l_spearman_coeff)):
        if l_p_value_1[i] < 0.05:
            f.write(str(l_spearman_coeff[i])+"\n")
        else:
           f.write(str(0) + "\n")

pos = [2 +i*0.5 for i in range(77)]
matrix_pat_array = np.array(matrix_pat)
matrix_score_array = np.array(matrix_score)
matrix_pat_array=np.transpose(matrix_pat_array,[1,0])
matrix_score_array=np.transpose(matrix_score_array,[1,0])
#
X = np.array([matrix_pat_array,matrix_score_array])
distances = distance_freq(pos)
threshold = 1.1
freq_adjacency = csr_matrix((distances <= threshold).astype(int))

thresh = 0.41 #0.41

T_obs_2, clusters_2, p_values_2, HO = mne.stats.permutation_cluster_test(X, threshold=thresh, stat_fun=perm_t_stat,
                                                                        n_permutations=1000, tail=0,
                                                                        adjacency=freq_adjacency, n_jobs=-1)

print(clusters_2)
print(p_values_2)
