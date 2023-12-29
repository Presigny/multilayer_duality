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

### Sripts

### Manual
Steps to reproduce the results of *Figure5a-b-c-d, Figure S4, Table S2, Table S3*
**Figure 5c**:




