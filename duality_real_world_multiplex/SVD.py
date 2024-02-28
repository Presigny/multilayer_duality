"""Author: Charley Presigny
mail: charley.presigny@inria.fr
licence: Creative Commons BY 4.0"""
import numpy as np
import pickle
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import scipy.sparse as sps


def to_supra_index(N,uplet):
    """Convert tensorial indices (layer,layer,node,node) into flattened matrix indices (row,col)
        Input:
        N -- number of nodes
        uplet -- tensor indices (layer,layer,node,node)
        Output:
        row, col -- matrix indices (row,col)"""
    k, l, i, j = uplet[0],uplet[1],uplet[2],uplet[3]
    row = k*N +i
    col = l*N +j
    return row,col

def to_tensor_index(N,uplet):
    """Convert flattened matrix indices (row,col) into tensorial indices (layer,layer,node,node)
        Input:
        N -- number of nodes
        uplet -- matrix indices (row,col)
        Output:
        row, col -- tensor indices (layer,layer,node,node) """
    row,col = uplet[0],uplet[1]
    k = row // N
    l = col // N
    i = row % N
    j = col % N
    return k,l,i,j

def compute_svd(A,N,M):
    """Compute the SVD decomposition of the contribution matrix
            Input:
            A -- sparse supra-adjency matrix
            N -- number of nodes
            M -- number of layers
            Output:
            C -- contribution matrix
            u -- left unitary matrix
            s -- singular value vectors
            v_t -- right unitary matrix """
    l_row, l_col = A.nonzero()
    C = np.zeros(shape=[N, M])
    for row, col in zip(l_row, l_col):
        k, l, i, j = to_tensor_index(N, [row, col])
        C[i, k] += A[row, col]
        C[j, l] += A[row, col]
    u, s, v_t = linalg.svd(C)
    return C,u,s,v_t

def to_Dark_side(s_A,N,M):
    """Generate the layerwise supra-adjacency matrix from the nodewise one
            Input:
            s_a -- sparse supra-adjacency matrix
            N -- number of nodes
            M -- number of layers
            Output:
            s_B -- sparse nodewise supra-adjacency matrix"""
    l_row, l_col = s_A.nonzero()
    s_B = sps.lil_matrix((M*N,M*N))
    num = 0
    for row, col in zip(l_row, l_col):
        num += 1
        k, l, i, j = to_tensor_index(N, [row, col])
        row_B,col_B = to_supra_index(M,[i,j,k,l]) # ''nodes become layers and layers become nodes''
        s_B[row_B,col_B] = s_A[row,col]
    return s_B

#####################################################################################################################
if __name__ == "__main__":
    ## Parameters
    path = "mplex_lil_sparse_random_multilayer_N_50_M_50.dat" #
    marker = "v"#"^"#"D"#"^"#"X"#"^"#"D"#"^"
    alpha_node = 0.5 #transparency of the nodewise vector lines
    thickness_node = 0.5 #thickness of the nodewise vector lines
    alpha_layer = 0.3 #transparency layerwise
    thickness_layer = 0.3 #thickness nodewise
    #percent = 75 # proportion of points to show in the inset
    multiplier = 10 # multiplier of e_tilde vector to reproduce it or not (see condition below)
    scale= "linear"

    ##Load data
    with open(path, 'rb') as f:
        parameter = pickle.load(f)  # we load the data file containing almost all the parameters already
    s_A = parameter['sparse']
    N = parameter['N']
    M = parameter['M']
    ## Run the SVD decomposition in the nodewise side
    C,u,s,v_t =  compute_svd(s_A,N,M)
    ## MAKE THE SINGULAR VECTOR POSITIVE
    u_1 = [abs(u[i][0]) for i in range(len(u))]
    u_2 = [abs(u[i][1]) for i in range(len(u))]
    v_t[0] = [abs(v_t[0][i]) for i in range(len(v_t[0]))]
    v_t[1] = [abs(v_t[1][i]) for i in range(len(v_t[1]))]
    s_real = np.zeros([N,M])
    for i in range(len(s)):
        s_real[i,i] = s[i]
    s_inv = linalg.pinv(s_real)
    s_2_inv_u1 = s_inv[:,0]
    s_2_inv_u2 = s_inv[:,1]
    ##Generate the vectors
    e_tilde_u1 = np.matmul(s_2_inv_u1,v_t) #the layer vectors projected on U1
    e_tilde_u2 =np.matmul(s_2_inv_u2,v_t) #the layer vector projected on U2
    n_tilde_u1 = np.zeros(N)
    n_tilde_u2 = np.zeros(N)
    for i in range(N):#=l_index:
        for k in range(M):
            n_tilde_u1[i] += C[i,k]*e_tilde_u1[k]
            n_tilde_u2[i] += C[i,k]*e_tilde_u2[k]

    n_tilde_u1_corrected = [n_tilde_u1[i] for i in range(len(n_tilde_u1)) if n_tilde_u1[i] != 0 or n_tilde_u2[i] != 0]
    n_tilde_u2_corrected = [n_tilde_u2[i] for i in range(len(n_tilde_u2)) if n_tilde_u2[i] != 0 or n_tilde_u1[i] != 0]
    ## Make the plot
    fig,ax1 = plt.subplots(figsize=(7,11))
    for i in range(M):
        if (multiplier*e_tilde_u1[i] > 1e-9 and multiplier*e_tilde_u1[i] > 1e-9) and e_tilde_u1[i] != 0 and e_tilde_u1[i] != 0: #reproduce or not the vector on the plot
            ax1.axline([0,0],[multiplier*e_tilde_u1[i],multiplier*e_tilde_u2[i]],linestyle='--',color="black",alpha=alpha_node,linewidth=thickness_node,zorder=3)
        ax1.scatter(n_tilde_u1_corrected, n_tilde_u2_corrected, color="#BDD7EE", alpha=1, edgecolors="black", s=300,marker=marker,zorder=5)

    plt.title("nodewise",fontsize=20)
    plt.xscale(scale)
    plt.yscale(scale)
    # Make the limits on x and y axes
    c = 1#0.12#0.8437294885207698#0.85#0.74#0.8437294885207698 #0.74
    d =1 #0.19#0.49502353123549714#2.70#0.84#0.49502353123549714 #0.84
    # ax1.set_xlim(-4e-2, 1.1 * c)
    # ax1.set_ylim(-9e-3, 1.1 * d)
    ax1.invert_xaxis()
    print(max(n_tilde_u1_corrected))
    print(max(n_tilde_u2_corrected))
    ax1.set_xlabel("U(1)",fontsize=17)
    ax1.set_ylabel("U(2)",rotation=0,labelpad=0,fontsize=17,y=1)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.savefig("Images/NODE_"+path+".png",format='png',dpi=300,transparent=True)
    plt.show()

    ##Convert the supra-adjacency matrix in the layerwise side
    s_B = to_Dark_side(s_A,N,M)
    ## Run the SVD decomposition in the layerwise side
    C,u,s,v_t =  compute_svd(s_B,M,N)
    ## MAKE THE SINGULAR VECTOR POSITIVE
    u_1 = [abs(u[i][0]) for i in range(len(u))]
    u_2 = [abs(u[i][1]) for i in range(len(u))]
    v_t[0] = [abs(v_t[0][i]) for i in range(len(v_t[0]))]
    v_t[1] = [abs(v_t[1][i]) for i in range(len(v_t[1]))]
    s_real = np.zeros([M,N])
    for i in range(len(s)):
        s_real[i,i] = s[i]
    s_inv = linalg.pinv(s_real)
    s_2_inv_u1 = s_inv[:,0]
    s_2_inv_u2 = s_inv[:,1]
    ##Generate the vectors
    e_tilde_u1 = np.matmul(s_2_inv_u1,v_t) #the layer vectors projected on U1
    e_tilde_u2 =np.matmul(s_2_inv_u2,v_t) #the layer vector projected on U2
    n_tilde_u1 = np.zeros(M)
    n_tilde_u2 = np.zeros(M)
    for i in range(M):#=l_index:
        for k in range(N):
            n_tilde_u1[i] += C[i,k]*e_tilde_u1[k]
            n_tilde_u2[i] += C[i,k]*e_tilde_u2[k]

    n_tilde_u1_corrected = [n_tilde_u1[i] for i in range(len(n_tilde_u1)) if n_tilde_u1[i] != 0 or n_tilde_u2[i] != 0]
    n_tilde_u2_corrected = [n_tilde_u2[i] for i in range(len(n_tilde_u2)) if n_tilde_u2[i] != 0 or n_tilde_u1[i] != 0]
    print(max(n_tilde_u1_corrected))
    print(max(n_tilde_u2_corrected))
    ##Make the plot
    fig,ax1 = plt.subplots(figsize=(7,11))
    for i in range(N): #in l_inde_layer
        if (multiplier*e_tilde_u1[i] > 1e-9 and multiplier*e_tilde_u1[i] > 1e-9) and e_tilde_u1[i] != 0 and e_tilde_u1[i] != 0:
            ax1.axline([0,0],[multiplier*e_tilde_u1[i],multiplier*e_tilde_u2[i]],linestyle='--',color="black",alpha=alpha_layer,linewidth=thickness_layer,zorder=3)
    ax1.scatter(n_tilde_u1_corrected,n_tilde_u2_corrected,color="#F8CBAD",alpha=1,edgecolors="black",s=300,marker=marker,zorder=5)
    plt.title("layerwise",fontsize=20)
    plt.xscale(scale)
    plt.yscale(scale)
    # ax1.set_xlim(-3e-2,1.1*c)
    # ax1.set_ylim(-1e-2, 1.1*d)
    ax1.set_xlabel("V(1)",fontsize=17)
    ax1.set_ylabel("V(2)",rotation=0,labelpad=0,fontsize=17,y=1.02)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.savefig("Images/LAYER_" + path+".png",format='png',dpi=300,transparent=True)
    plt.show()

