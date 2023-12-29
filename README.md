# About

This repository is meant to reproduce the results and figures from the
**Paper:** [1] Topological duality of multilayer networks, C. Presigny,MC. Corsi and F. De Vico Fallani,*arXiv*, 2023, https://arxiv.org/abs/2306.12136

**Version**:
- All python codes were tested with Python 3.9.0 or Jupyter 6.2.0 on macOS Monterrey 12.7.1
- *script_explain_distances_Fig3.m* was tested with MATLAB R2023a on macOS Monterrey 12.7.1

**Requirements:**
Required python packages to excecute the scripts:
- kneed (tested on version 0.8.3)
- matplotlib (tested on version 3.5.2)
- mne (tested on version 0.23.4)
- numpy (tested on version 1.23.0)
- pandas (tested on version 1.4.3)
- scipy (tested on version 1.9.1)
- scikit-learn (tested on version 1.0.1)
- statsmodels (tested on version 0.13.2)

**Directories:**
The scripts are organized into 3 directories associated to results in the paper:
- node_layer_duality: scripts to reproduce results of *Figure 2b-c,Figure 3, Figure S1* and multilayer stochastic rewiring algorithm
- duality_real_world_multiplex: scripts to process and analyze the real-world data for results corresponding to *Figure 4a-b, Table S1,Figure S3*
- duality_brain_networks: scripts to process and analyze the multifrequency brain networks for results corresponding to *Figure5a-b-c-d, Figure S4, Table S2, Table S3*

## node_layer_duality

## duality_real_world_multiplex

### Data

### Sripts

### Manual

## duality_brain_networks

This sections is dedicated to how to download, process and analyze the data that support the results corresponding to *Figure5a-b-c-d, Figure S4, Table S2, Table S3*
**Directories**:
- preprocessing: contains the script to process the supra-adjacency matrice before analyzing them *per se*
- preprocessing/database: contains the supra-adjacency matrices and the results of the analysis
- preprocessing/database: contains the reduced version supra-adajcency matrices (in terms of layers) and the results of the analysis
- metadata: contains metadata necessary for the analysis (distance between brain regions, name of brain regions, cognitive score of patients)


### Data

### Scripts

### Manual
Step-by-step to reproduce the results of *Figure5a-b-c-d, Figure S4, Table S2, Table S3*
**Figure 5c**:
1. Symmetrize the supra-adjency matrices, compute node and layer multistrength centrality for patients and healthy subjects separately 
   1. preprocessing/matrix_symmetrization.py: load_path = "database/individual_matrix/";  save_path = "database/sym_individual_matrix/"
   2. preprocessing/compute_multistrength.py: load_path = "database/sym_individual_matrix/"; dir_path = "database/sym_multistrength_PAT"; l_load_path = l_load_pat
   3. preprocessing/compute_multistrength.py: load_path = "database/sym_individual_matrix/"; dir_path = "database/sym_multistrength_SUJ"; l_load_path = l_load_suj
2. Generate supra-adjacency matrix averaged over healthy subjects, compute its node and layer multistrength centrality
   1. preprocessing/average_over_subjects.py: load_path = "database/sym_individual_matrix/"; save_path = "database/sym_average_SUJ_matrix/sym_average_SUJ.dat"
   2. preprocessing/compute_multistrength.py: load_path = "database/sym_individual_matrix/"; dir_path = "database/sym_multistrength_avgSUJ"; l_load_path = l_load_suj
3. Compute the multistrength-based euclidean distance between patients and average healthy
   1. euclidean_distance.py: main_path = "preprocessing/database/"patient_path = main_path+"sym_multistrength_PAT/"; load_path_avg = main_path+"/sym_multistrength_avgSUJ/sym_average_SUJ.dat_multistrength"; writing_path = main_path+"70_euclidean_distance_to_avg_SUJ_sym.csv"
4. Compute the reduced matrix and related computations
   1. preprocessing/reduce_matrix.py: Select the desired *l_sup_layer* list; dir_path = "database_reduced/%s_sym_individual_matrix/"%(number_layers)
   2. reproduce step 1,2,3 with the reduced matrices
5. Plot the Figure 5c
   1. plotter_distance_func_M.py

**Figure 5a-b**:
1. Brain_analysis_5ab.ipynb: compute the Wilcoxon rank-sum statistic between the node/layer multistrength distributions of patients and healthy subjects

**Figure 5d, Table S3**:
1. correlation_mmse_freq.py: Make plot MMSE score in function of multistrength for each considered frequency, compute their Spearman coefficient and related p values, correct the p values for multiple comparisons (Benjamini-Hochberg), correct p values with cluster-based permutation test
   1. save_path = "preprocessing/database/plot_freq_correlation_MMSE/": where the MMSE in fucntion of multistrength plot are saved
   2. l_spearman_coeff;l_p_value: where the Spearman coefficients and p_values of all MMSE, frequency multistrength correlation are stored
   3. writing_path_FDR = "preprocessing/database/FDR_0.05_freq_correlation_MMSE.csv": where the correction for multiple comparison is saved
   4. clusters_2;p_values_2: where the cluster and corrected p_values with cluster-based permutation test are stored

**Table S2**:
1. correlation_mmse_ROI.py: Make plot MMSE score in function of multistrength for each considered region, compute their Spearman coefficient and related p values, correct the p values for multiple comparisons (Benjamini-Hochberg), correct p values with cluster-based permutation test
   1. save_path = "preprocessing/database/plot_ROI_correlation_MMSE/": where the MMSE in function of multistrength plot are saved
   2. l_spearman_coeff;l_p_value: where the Spearman coefficients and p_values of all MMSE, frequency multistrength correlation are stored
   3. writing_path_FDR = "preprocessing/database/FDR_0.05_ROI_correlation_MMSE.csv": where the correction for multiple comparison is saved
   4. clusters_2;p_values_2: where the cluster and corrected p_values with cluster-based permutation test are stored

**Table S4**:
1. Compute and plot the absolute difference between the intralayer, interlayer and replica weights of patients matrices and the supra-adjacency matrix averaged over healthy subjects
  1. check_diff_weights_S4.py: load_path = "preprocessing/database_reduced/%s_sym_individual_matrix_PAT/"%(number_of_layer); load_HC_average = "preprocessing/database_reduced/%s_sym_average_SUJ_matrix/%s_sym_average_SUJ.dat"%(number_of_layer,number_of_layer)

