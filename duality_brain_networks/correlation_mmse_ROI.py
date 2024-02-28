"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import mne
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests


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
        a = dataframe_2[dataframe_2['Code_IRM'].str.contains(name[0])] #find the patient in dataframe_2
        if len(a) != 0:
            dataframe_1.at[ind, score] = int(a[score])
        else:
            dataframe_1.at[ind, score] = -1 #put -1 if no score is associated to a given patient
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
    l_tvalue = np.zeros([len(m_pat[0])])
    l_pvalue = np.zeros([len(m_pat[0])])
    for i in range(len(m_pat[0])):
        t_values, pv = spearmanr(m_pat[:,i],l_score)
        l_tvalue[i] = t_values
        l_pvalue[i] = pv
    return l_tvalue

##BEGINNING OF THE SCRIPT
l_p_value = [] #list where the p_values of each correaltion are stored
l_spearman_coeff = [] #list where the spearman coeff is stored
matrix_pat = []
matrix_score = []
path_node_label = "metadata/name_of_ROIS.txt"
data_anat_labels = np.genfromtxt(path_node_label, dtype='str')[0:,0] #get the names of ROIs
l_ROI_index = range(70) #70 is the number of ROI cosniderd in the study

for ROI_index in l_ROI_index: #Do the process for each ROI
    layer_index = 66
    writing_path = "preprocessing/database/ROI_Activity_PAT.csv"
    restore_path = "preprocessing/database/sym_multistrength_PAT/"
    save_path = "preprocessing/database/plot_ROI_correlation_MMSE/"
    l_restore_path = os.listdir(restore_path)

    l_name = [l_restore_path[i][26:-18] for i in range(len(l_restore_path))]
    with open(writing_path, 'w') as f:  # save name,N,M,std_X,std_Y,density
            string_csv = "name,K_X_"+str(ROI_index)+",K_Y_"+str(layer_index)+"\n"
            f.write(string_csv)
    for patient_path in l_restore_path:
        with open(restore_path + patient_path, 'rb') as f: #open file 1
            dico = pickle.load(f)
        multistrength = dico['multistrength'][ROI_index]
        multistrength_DS = dico['multistrength_DS'][layer_index]
        with open(writing_path, 'a') as f:  # save name,N,M,std_X,std_Y,density
            string_csv = str(patient_path[26:-18]) + "," + str(multistrength) + "," + str(multistrength_DS) + "\n"
            f.write(string_csv)
    # Load the excel file containing the MMSE score of patients
    df_score = pd.read_excel('metadata/PAT_metadata.xlsx')
    df_activity = pd.read_csv(writing_path, encoding='unicode_escape', index_col='name')
    print(df_score.columns)
    df_score = df_score[['groupe','MMSE','Code_IRM','âge','Sexe']]
    for i in df_score.index: #change the name
        df_score['Code_IRM'] = df_score['Code_IRM'].replace(df_score['Code_IRM'][i], df_score['Code_IRM'][i][17:])
    # Separate patients from subjects
    df_score_pat = df_score[df_score['Code_IRM'].str.contains('AT')]
    df_score_suj = df_score[df_score['Code_IRM'].str.contains('SUJ|Sujet')]

    df_score_pat = clean_scores(df_score_pat)
    score = 'MMSE'
    df_activity = score_in_df(df_activity, df_score_pat, 'MMSE')
    l_KX,l_KY,l_score = recalibrate(df_activity,score)
    matrix_pat.append(l_KX)
    matrix_score.append(l_score)
    #l_dX,l_dY,l_score,l_outlier_X,l_outlier_Y = remove_outliers(l_dX,l_dY,l_score)

    fontsize=18
    # save the correlation between MMSE score and the total activity
    correlation = spearmanr(l_KX,l_score)
    coeff = round(correlation[0], 2)
    pvalue = round(correlation[1],3)
    ##PLOT
    if (ROI_index%2) == 0: #select if it is the left or right regions
        u ="L "
    else:
        u="R "
    first = r"$S$ = %s , $p$ = %s  "%(coeff,pvalue)
    second=r"ROI = "+u+data_anat_labels[ROI_index]
    label = first+second
    plt.scatter(l_KX,l_score,color='#BDD7EE',edgecolor='black',s=100,alpha=1,label=label,zorder=5)
    plt.xlabel(r"$k_\mathcal{X}^{AD}$",fontsize=fontsize,labelpad=-3)
    plt.ylabel("MMSE",rotation=90,fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax1 = plt.gca()
    ax1.set_aspect(1/ax1.get_data_ratio(),adjustable='box')
    plt.savefig(save_path + "ROI_" + str(ROI_index) + "_article.png", format='png')
    #plt.show()
    plt.close()
    l_spearman_coeff.append(correlation[0])
    l_p_value.append(correlation[1])


## Do the FDR correction on the list of p value for each frequency to correct for multiple comparisons

a,b,c,d = multipletests(l_p_value, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)


##Do the cluster based permutation test to correct for multiple comparisons
matrix_pat_array = np.array(matrix_pat)
matrix_score_array = np.array(matrix_score)
matrix_pat_array = np.transpose(matrix_pat_array, [1, 0])
matrix_score_array = np.transpose(matrix_score_array, [1, 0])
#
X = np.array([matrix_pat_array, matrix_score_array])
distances = np.load("metadata/euclid_coo.npy") #distance matrices between ROIS
threshold = 1
ROI_adjacency = csr_matrix((distances <= threshold).astype(int))

thresh = 0.41 #0.41

T_obs_2, clusters_2, p_values_2, HO = mne.stats.permutation_cluster_test(X, threshold=thresh, stat_fun=perm_t_stat,
                                                                        n_permutations=1000, tail=0,
                                                                        adjacency=ROI_adjacency, n_jobs=-1)
## The obtained cluster and their related corrected p_values
print(clusters_2)
print(p_values_2)
