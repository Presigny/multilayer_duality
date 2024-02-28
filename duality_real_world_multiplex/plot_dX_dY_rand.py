"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import pandas as pd
import os
import sklearn.cluster as sk
import itertools
import matplotlib.pyplot as plt
from kneed import KneeLocator
from sklearn.metrics import silhouette_score
from scipy.stats import pearsonr

def plotter_cluster_label(l_x,l_y,label_x,label_y,cluster_label,scale='linear'):
    """ Plot dY/dY_rand in function of dX/dX_rand for several group of empirical clusters.
    Colors reflects their clustering assignement according to a clustering algorithm
            Input:
            l_x -- log(dX/dX_rand) for each empirical network in each group
            l_y -- log(dY/dY_rand) for each empirical network in each group
            label_x -- name of x-axis
            label_y -- name of the y-axis
            cluster_label -- cluster assignement of each empirical network
            scale -- scale of the plot
            Output:
            Show and save the plot
            """
    fig,ax = plt.subplots(figsize=(20.5, 22))
    plt.rcParams.update({'font.size': 38})
    l_label = ["TwitterEvent", "PierreAuger", "Arxiv", "Genetic", "HumanMicrobiome", "FAOTrade",
                    "UgandaVillage", "GermanTransport", "EUAir", "C.Elegans"]  # "spatial_socio_eco"

    for i in range(len(l_y)):
        for j in range(len(l_y[i])):
            if j == 0:
                plt.scatter(l_x[i][j], l_y[i][j], label=l_label[i], marker=mStyles[i],
                            color=l_color[cluster_label[i][j]], s=400, alpha=1, edgecolors='black')
            else:
                plt.scatter(l_x[i][j], l_y[i][j], marker=mStyles[i],
                            color=l_color[cluster_label[i][j]], s=400, alpha=1, edgecolors='black')

    plt.grid()
    plt.legend(fontsize=32)
    plt.yticks(fontsize=32)
    plt.xticks(fontsize = 32)
    plt.ylabel(label_y, size=48,rotation=0,labelpad=20,y=0.95)
    plt.xlabel(label_x, size=48,x=1,labelpad=-20,y=0.2)
    plt.xscale(scale)
    plt.yscale(scale)
    ax.spines[['right', 'top']].set_visible(False)
    plt.savefig("main_fig_clustering_mplex_no_line.png",format='png',dpi=300,transparent=True)
    print("OK")
    plt.show()
    plt.close()
    return 0

if __name__ == "__main__":
    ##Select the normalized distances associated with each empirical datasets
    os.chdir("normalized_distance")
    file_to_read = ["twitter_event","PierreAuger","Arxiv","genetic","HumanMicrobiome","FAO_trade","Uganda_villages","german_transport","EUAir","C_Elegans"] #"spatial_socio_eco"
    #
    l_x,l_y,l_avg_Kx,l_avg_Ky,l_density,l_N,l_M,l_name = [],[],[],[],[],[],[],[]
    df = pd.DataFrame()
    l_label = []

    ## Load all the data into the required lists
    for file in file_to_read:
        df_temp = pd.read_csv(file)
        df_temp = df_temp.sort_values(by=['std_y'], ascending=True)  # sort the value from the normalized layerwise distance
        l_y.append(df_temp["std_y"].values.tolist()) #store the normalized layerwise distance (d_Y/d_Y_rand)
        l_x.append(df_temp["std_x"].values.tolist()) #store the normalized nodewise distance (d_X/d_X_rand)
        l_avg_Ky.append(df_temp["K_y"].values.tolist())
        l_avg_Kx.append(df_temp["K_x"].values.tolist())
        l_density.append(df_temp["density"].values.tolist()) #store the density of each empirical network
        l_N.append(df_temp["N"].values.tolist()) #store the number of nodes
        l_M.append(df_temp["M"].values.tolist()) #store the number of layers
        l_label.append(file) #store the name of the group of empirical network
        l_name .append(df_temp["name"].values.tolist()) #store the names of each empirical network of a group
        df = df.append(df_temp)
    ## Parameters
    mStyles = ["D","v","^","P",">","s","p","D","d","X","H","+","x","h","D","d","|","_",0,1,2,3,4,5,6,7,8,9,10,11] #style of markers for the plot
    spatial_color = '#7A499F' #color of the
    n_cluster = 2 #choose the number of clusters to compute

    ## Formatting and creating variables for the script

    l_x_rand_log = [np.log(l_x[i]) for i in range(len(l_x))] #to have data in log scale
    l_y_rand_log = [np.log(l_y[i]) for i in range(len(l_y))]
    label_cluster =  [] #cluster assignement of each empirical network
    for i in range(len(l_x)):
        label_cluster.append([0]*len(l_x[i]))
    # #
    l_X = [] #make couple (dX,dY) for every empirical network
    l_x_1D = []
    l_y_1D = []
    for i in range(len(l_x)):
        for x,y in zip(l_x_rand_log[i],l_y_rand_log[i]):
            l_X.append([x,y])
            l_x_1D.append(x)
            l_y_1D.append(y)

    ## Clustering using Kmeans clustering algorithm,
    kmeans = sk.KMeans(n_clusters=n_cluster,tol=1e-8,max_iter=1000).fit(l_X)
    l_result = kmeans.labels_ #the assignement of each empirical network
    inertia = kmeans.inertia_ #variable to evaluate the quality of the cluster
    h = iter(l_result) # associate the cluster assignement of each network as
    for i in range(len(l_x)):
        for j in range(len(l_x[i])):
            label_cluster[i][j] = next(h)
    ## Clustering using Agglomerative clustering algorithm,
    # kmeans = sk.AgglomerativeClustering(n_clusters=n_cluster).fit(l_X)
    # l_result = kmeans.labels_
    # #inertia = kmeans.inertia_
    # h = iter(l_result)
    # for i in range(len(l_x)):
    #     for j in range(len(l_x[i])):
    #         label_cluster[i][j] = next(h)

    l_color = [spatial_color,'#96BF0D', '#25fde9','b','g','r','c','m','y',"k"]
    l_color = ['#96BF0D',spatial_color,'#25fde9'] #list of color for each group of empirical networks
    plotter_cluster_label(l_x_rand_log,l_y_rand_log,r"$\frac{d_X}{d_X^{ER}}$",r"$\frac{d_Y}{d_Y^{ER}}$",label_cluster,scale="linear")

    ## Elbow test for clustering quality
    sse = []
    for k in range(1, 11): #use the kmean clustering from k=1 cluster to k=11 clusters
         kmeans = sk.KMeans(n_clusters=k,tol=1e-8,max_iter=1000).fit(l_X)
         sse.append(kmeans.inertia_) #save the sum of squared error
    # Plot the of the elbow method
    fig, ax = plt.subplots(figsize=(10,10))
    plt.style.use("fivethirtyeight")
    plt.plot(range(1, 11), sse)
    plt.xticks(range(1, 11))
    plt.xlabel("Number of Clusters")
    plt.ylabel("SSE")
    ax.set_aspect(1. / ax.get_data_ratio(), adjustable='box')
    plt.savefig("kmeans_elbow_score.png",format="png",transparent=True)
    plt.show()
    plt.close()

    kl = KneeLocator(range(1, 11), sse, curve="convex", direction="decreasing") #locate the knee point on the figure
    print(kl.elbow)

    ##Silhouette score and plot for Kmeans clustering

    # A list holds the silhouette coefficients for each k
    silhouette_coefficients = []
    # Notice you start at 2 clusters for silhouette coefficient
    for k in range(2, 11):
         kmeans = sk.KMeans(n_clusters=k,tol=1e-8,max_iter=1000).fit(l_X)
         score = silhouette_score(l_X, kmeans.labels_)
         silhouette_coefficients.append(score)

    fig, ax = plt.subplots(figsize=(10,10))
    plt.style.use("fivethirtyeight")
    plt.plot(range(2, 11), silhouette_coefficients)
    plt.xticks(range(2, 11))
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Coefficient")
    ax.set_aspect(1. / ax.get_data_ratio(), adjustable='box')
    plt.savefig("kmeans_silhouette_score.png",format="png",transparent=True)
    plt.show()
    #
    ##Silhouette score and plot for agglomerative clustering

    # A list holds the silhouette coefficients for each k
    silhouette_coefficients = []
    # Notice you start at 2 clusters for silhouette coefficient
    for k in range(2, 11):
         kmeans = sk.AgglomerativeClustering(n_clusters=k).fit(l_X)
         score = silhouette_score(l_X, kmeans.labels_)
         silhouette_coefficients.append(score)

    fig, ax = plt.subplots(figsize=(10,10))
    plt.style.use("fivethirtyeight")
    plt.plot(range(2, 11), silhouette_coefficients)
    plt.xticks(range(2, 11))
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Coefficient")
    ax.set_aspect(1. / ax.get_data_ratio(), adjustable='box')
    plt.savefig("agglo_silhouette_score.png",format="png",transparent=True)
    plt.show()


## Correlation dY/DX in the German transport dataset

    print(pearsonr(l_x_rand_log[7], l_y_rand_log[7])) #Pearson correlation
    l_x_german = np.array(l_x_rand_log[7]).reshape(-1, 1)
    l_y_german = np.array(l_y_rand_log[7]).reshape(-1, 1)

    import sklearn.linear_model as sklin

    reg = sklin.LinearRegression().fit(l_x_german, l_y_german) #Linear regression
    l_predicted_curve = reg.predict(l_x_german)
    l_R_square = reg.score(l_x_german, l_y_german)
    l_slope = reg.coef_[0]