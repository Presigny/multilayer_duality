# About

This repository is meant to reproduce the results and figures from the 
**Paper:** \ [[1]](#1) Topological duality of multilayer networks, C. Presigny, MC. Corsi and F. De Vico Fallani,*arXiv*, 2023, https://arxiv.org/abs/2306.12136

**Version**:
- All python codes were tested with Python 3.9.0 or Jupyter 6.2.0 on macOS Monterrey 12.7.1
- *script_explain_distances_Fig3.m* was tested with MATLAB R2023a on macOS Monterrey 12.7.1

**Requirements:** \
Required python packages to excecute the scripts:
- kneed (tested on version 0.8.3)
- matplotlib (tested on version 3.5.2)
- mne (tested on version 0.23.4)
- numpy (tested on version 1.23.0)
- pandas (tested on version 1.4.3)
- scipy (tested on version 1.9.1)
- scikit-learn (tested on version 1.0.1)
- statsmodels (tested on version 0.13.2)

**Directories:** \
The scripts are organized into 3 directories associated to results in the paper:
- node_layer_duality: scripts to reproduce results of *Figure 2b-c,Figure 3, Figure S1* and multilayer stochastic rewiring algorithm
- duality_real_world_multiplex: scripts to process and analyze the real-world data for results corresponding to *Figure 4a-b, Table S1,Figure S3*
- duality_brain_networks: scripts to process and analyze the multifrequency brain networks for results corresponding to *Figure5a-b-c-d, Figure S4, Table S2, Table S3*

## node_layer_duality

This sections is dedicated to how to rewire numerically a multilayer network according to the stochastic rewiring model presented in 1, how to generate random multilayer and multiplex networks, how to process and analyze the data that support the results corresponding to *Figure 2b-c, Figure 3, Figure S1*\
**Directories**: \
- mock_data: contains the mock_data to use the scripts to generate *Figure 2b-c, Figure 3, Figure S1*
- mock_data/distance_diagramme: contains the data to generate *Figure 2c*
- mock_data/full_N_200_M_200: contains the random multilayer networks (number of nodes N=200, number of layers M=200) to generate *Figure 2b*
- mock_data/mplex_N_200_M_200: contains random multiplex networks (number of nodes N=200, number of layers M=200)
- mock_data/sparse_rewiring_data: contains the data to generate *Figure S1*

### Scripts
- full_multilayer_fig2b.py: uses the data in mock_data/full_N_200_M_200 to generate the node and layer distances in function of the probability parameters and proportion of rewiring R.
- analysis_of_rewiring_figS1.py: plot the data contained in mock_data/sparse_rewiring_data to obtain *Figure S1* (comparison between numerical and theoretical distances between expected and initial node/layer multidegree centrality distances
- distance_diagramme.py: plot the data contained in mock_data/distance_diagramme to obtain *Figure 2c* (theoretical distances in function of the probability parameters p_node and p_layer)
- full_model_fig2c.py: uses the data in mock_data/full_N_200_M_200 to generate the theoretical node and layer distances in function of the probability parameters (p_node,p_layer) for multilayer networks
- mplex_model_fig2c.py: uses the data in mock_data/mplex_N_200_M_200 to generate the theoretical node and layer distances in function of the probability parameters (p_node,p_layer) for multiplex networks
- plotter_d_Y_in_function_of_d_X.py: plot the data generated in full_multilayer_fig2b.py to obtain *Figure 2b* (comparison between numerical and theoretical distances between expected and initial node/layer multidegree centrality distances
- script_explain_distances_Fig3.m: plot the the distances of random multiplex and random multilayer networks to obtain *Figure 3*
- sparse_random_generator_git.py: generate random multiplex/multilayer networks (number of nodes/layers and probability parameters can be set)
- sparse_rewiring.py: rewire numerically a multilayer/multiplex network according to the probability parameters p_node, p_layer, p_tel. Return the rewired supra adjancency matrix and the related node/layer multidegree sequence. Works only for binary networks.

### Manual
Step-by-step to reproduce the data in mock_data (not mentionned in the Scripts section) \
**mock_data/full_N_200_M_200** \
sparse_random_generator.py: n = 100; N = 200; M = 200;density = 0.005; interlayer = True; multiplex_random = True; multiplex = False;
path = "mock_data/full_N_200_M_200/mock_random_multilayer_network_"+str(n) \
**mock_data/mplex_N_200_M_200** \
sparse_random_generator.py: n = 100; N = 200; M = 200;density = 0.05; interlayer = False; multiplex_random = False; multiplex = True;path = "mock_data/mplex_N_200_M_200/mock_random_multiplex_network_"+str(n) \
**mock_data/sparse_rewiring_data** \
sparse_rewiring.py: uncomment idx, save_path_original and save_path_DS (line 207-209);uncomment lines 233-238; put idx = 0 when event_proba = [1, 0, 0]; put idx = 64 when event_proba = [0, 1, 0]; put idx = 74 when event_proba = [0, 0, 1]; change R to have a different proportion of rewiring

## duality_real_world_multiplex
This sections is dedicated to how to download, process and analyze the data that support the results corresponding to *Figure 4a-b, Table S1,Figure S3* \
**Directories**: \
- preprocessing: contains the script to generate the supra-adjacency matrices of the multiplex networks from raw data
- preprocessing/database: contains the preprocessed multiplexes
- normalized_distance: contains the normalized distance between multiplexes and rewired version of them.Needed to produce Fiure 4a


### Data
Details on the data used in [[1]](#1) can be retrieved in the associated Supplementary material. For simplicity, data are preprocessed already in the repository even if the code to do so is available as well.
- Twitter events, PierreAuger, Arxiv, EuAir, Genetic, C.Elegans, FAO trade, HumanMicrobiome: datasets are available [here](https://manliodedomenico.com/data.php) [[3-11]](#3)
- German transport: Refer to [Urban-multiplex-networks repository](https://github.com/KBergermann/Urban-multiplex-networks) to find the methods to generate and download the transport multiplex networks from the GTFS data [[12]](#12)
- Uganda villages: dataset is available [here](https://doi.org/10.17863/CAM.15616) [[12]](#12)

### Scripts
- preprocessing/pre_processing_de_domenico.py: Generate the supra_adjacency matrices and related multidegree centrality sequences of the considered datasets of this [database](https://manliodedomenico.com/data.php) (see section Data) .Make sure to input the right number of nodes and layers before lauching the script
- preprocessing/pre_processing_german_transport.py: Generate the supra_adjacency matrices and related multidegree centrality sequences of the dataset available [here](https://github.com/KBergermann/Urban-multiplex-networks).
- preprocessing/pre_processing_villages.py: Generate the supra_adjacency matrices and related multidegree centrality sequences of the Uganda villages dataset (see section Data) .Make sure to input the right number of nodes and layers before lauching the script
- SVD.py: uses the data of a multilayer network to produce the singular vectors to project its contribution matrix into a 2D plane [[2]](#2). Make the plot of *Figure 4b* . Adjust c,d to set limits of the x and y axes
- compute_std_tabS1.py: compute the standard deviation of the node and layer multidegree centeality sequences and store them in a csv file. It allow to produce *Table S1*
- generate_fitting_random_multidegree.py: generate the multidegree sequences of random multiplex networks with the same density as the input multiplex network. Compute the distance between the expected multidegree centrality sequences of the input multiplex network from the rewiring and the ones of the input multiplex networks. Compute those distances for the generated sequences of random multiplexes. Save in a file the ratio between such distance of the input multiplex network and the average of the distances of the random multiplex networks *i.e.* the normalized distances.
- plot_dX_dY_rand.py: uses the data from generate_fitting_random_multidegree.py to plot the layer normalized distances in function of the node normalized distances. Cluster the multiplex network according to either Kmeans or agglomerative clustering. Make the elbow test and silhouette score of Kmeans clustering and plot them to obtain *Figure S3*. Plot the silhouette score for agglomerative clustering.


## duality_brain_networks

This sections is dedicated to how to download, process and analyze the data that support the results corresponding to *Figure5a-b-c-d, Figure S4, Table S2, Table S3* \
**Directories**: \
- preprocessing: contains the script to process the supra-adjacency matrice before analyzing them *per se*
- preprocessing/database: contains the supra-adjacency matrices and the results of the analysis
- preprocessing/database: contains the reduced version supra-adajcency matrices (in terms of layers) and the results of the analysis
- metadata: contains metadata necessary for the analysis (distance between brain regions, name of brain regions, cognitive score of patients)


### Data

### Scripts

### Manual
Step-by-step to reproduce the results of *Figure5a-b-c-d, Figure S4, Table S2, Table S3* \
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


## References
<a id="1">[1]</a>
Presigny C.,Corsi MC. and De Vico Fallani F. (2023)
Topological duality of multilayer networks. ,*arXiv*, 
arXiv, https://arxiv.org/abs/2306.12136

<a id="2">[2]</a> 
Dijkstra, E. W. (1968). 
Go to statement considered harmful. 
Communications of the ACM, 11(3), 147-148.
Arenas A., Borge-Holthoefer J., Gomez S., Zamora-Lopez G. (2010)
Optimal map of the modular structure of complex networks. 
New Journal of Physics, 12(5)

<a id="3">[3]</a> 
Omodei E., De Domenico M., Arenas A. (2015). 
Characterizing interactions in online social networks during exceptional events. 
Frontiers in Physics, 3

<a id="4">[4]</a> 
De Domenico M., Altmann E.G. (2020). 
Unraveling the Origin of Social Bursts in Collective Attention. Scientific Reports. 
Scientific Reports, 10(1)

<a id="5">[5]</a> 
De Domenico M., Lancichinetti A., Arenas A, Rosvall M. (2015).
Identifying Modular Flows on Multilayer Networks Reveals Highly Overlapping Organization in Interconnected Systems.
Physical Review X, 5(1)

<a id="6">[6]</a> 
Cardillo A., Gomez-Gardenes J., Zanin M., Romance M., Papo D., Pozo F., et al. (2013)
Emergence of network features from multiplexity. 
Scientific Reports, 3(1)

<a id="7">[7]</a> 
De Domenico M., Nicosia V., Arenas A., Latora V. (2015)
Structural reducibility of multilayer networks. 
Nature Communications, 6(1)

<a id="8">[8]</a> 
Chen BL., Hall DH., Chklovskii DB. (2006)
Wiring optimization can relate neuronal structure and function. 
Proceedings of the National Academy of Sciences, 103(12)

<a id="9">[9]</a> 
De Domenico M., Porter MA., Arenas A. (2015)
MuxViz: a tool for multilayer analy- sis and visualization of networks. 
Journal of Complex Networks, 3(2)

<a id="10">[10]</a> 
Ding T., Schloss PD. (2014)
Dynamics and associations of microbial community types across the human body.
Nature, 509(7500):357–360

<a id="11">[11]</a>
De Domenico M., Biamonte J. (2016)
Spectral Entropies as Information-Theoretic Tools for Complex Network Comparison.
Physical Review X., 6(4)

<a id="12">[12]</a>
Bergermann K., Stoll M. (2021)
Orientations and matrix function-based centralities in multiplex network analysis of urban public transport.
Applied Network Science,6(1):90

<a id="13">[13]</a> 
Chami GF., Ahnert SE., Kabatereine NB., Tukahebwa EM. (2017)
Social network fragmentation and community health.
Proceedings of the National Academy of Sciences, 114(36)

