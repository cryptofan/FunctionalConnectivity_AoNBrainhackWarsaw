# This file is intended to work as a library.
# you can
# import fc_methods as fc
# and then e.g.:
# fc.manhattan_distance(x,y)
# etc.
# Happy coding! :)


# Import dependencies
from scipy import special, spatial#, stats, linalg
from scipy.integrate import odeint
import numpy as np
from scipy import stats, linalg

# Define methods
def pearson_corr(matrix):
    """
    Function returns a matrix with Pearson correlation coefficients
    between pairs of variables (columns of the matrix) along with p-values.
    
    Parameters
    ----------
    matrix : array-like, shape (n, m)
             Array with each column being a different variable.
             
    Returns
    -------
    corr : array-like, shape (m, m)
           corr[i, j] contains the partial correlation of matrix[:, i] and matrix[:, j].
    prob : array-like, shape (m, m)
           prob[i, j] contains the p-value of a coefficient in corr[i, j].
    """
    (n, m) = matrix.shape

    DO = matrix - (np.sum(matrix, 0) / np.double(n))
    # note that mean row will be applyed row-wise to original matrices
    
    cov = np.einsum("nt,nm->tm", DO, DO)

    varO = np.sum(DO ** 2, 0)
    tmp = np.outer(varO, varO)
    
    corr = cov / np.sqrt(tmp)
    
    df = n-2
    
    diag = (np.diag(np.diag(corr)))
    corr -= diag
    t_squared = corr*corr*(df / ((1.0 - corr) * (1.0 + corr)))
    prob = special.betainc(0.5*df, 0.5, df / (df + t_squared))
    np.fill_diagonal(corr, 1)
    np.fill_diagonal(prob, 0)
    
    return corr, prob


def partial_corr(matrix):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in a matrix,
    controlling for the remaining variables in that matrix.
    
    Parameters
    ----------
    matrix : array-like, shape (n, m)
             Array with the different variables. Each column of the matrix is taken as a variable.
    
    Returns
    -------
    partial : array-like, shape (m, m)
              partial[i, j] contains the partial correlation of matrix[:, i] and matrix[:, j],
              controlling for the remaining variables in the matrix.
    prob : array-like, shape (m, m)
           prob[i, j] contains the p-value of a coefficient in corr[i, j].
    """

    n = matrix.shape[0]
    m = matrix.shape[1]
    ic = -np.linalg.pinv(np.cov(matrix, rowvar=0))
    diag1 = np.tile(np.sqrt(np.abs(np.diag(ic))),[m,1]).T
    diag2 = np.tile(np.sqrt(np.abs(np.diag(ic))),[m,1])
    partial = ((ic/diag1)/diag2)+2*np.eye(m)
    
    if n > m:
        df = n-m
    
        diag = (np.diag(np.diag(partial)))
        partial -= diag
        t_squared = partial*partial*(df / ((1.0 - partial) * (1.0 + partial)))
        prob = special.betainc(0.5*df, 0.5, df / (df + t_squared))
        np.fill_diagonal(partial, 1)
        np.fill_diagonal(prob, 0)
        return partial, prob
    else:
        return partial
    

def partial_partial_corr(matrix, alpha):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in a matrix,
    controlling for the remaining variables in that matrix, with additional option to partially partial out the variables.
    
    Parameters
    ----------
    matrix : array-like, shape (n, m)
             Array with the different variables. Each column of the matrix is taken as a variable.
    alpha :  parameter controlling the extent to which remaining variables are partialled out (0 - not at all, 1 - fully)
    
    Returns
    -------
    partial : array-like, shape (m, m)
              partial[i, j] contains the partial correlation of matrix[:, i] and matrix[:, j],
              controlling for the remaining variables in the matrix.
    prob : array-like, shape (m, m)
           prob[i, j] contains the p-value of a coefficient in corr[i, j].
    """
 
    n = matrix.shape[0]
    m = matrix.shape[1]
    partial = np.zeros((m, m), dtype=np.float)
    for i in range(m):
        partial[i, i] = 1
        for j in range(i+1, m):
            idx = np.ones(m, dtype=np.bool)
            idx[i] = False
            idx[j] = False
            beta_i = linalg.lstsq(matrix[:, idx], matrix[:, j])[0]
            beta_j = linalg.lstsq(matrix[:, idx], matrix[:, i])[0]

            res_j = matrix[:, j] - alpha*matrix[:, idx].dot(beta_i)
            res_i = matrix[:, i] - alpha*matrix[:, idx].dot(beta_j)

            corr = stats.pearsonr(res_i, res_j)[0]
            partial[i, j] = corr
            partial[j, i] = corr

    if n > m:
        df = n-m
    
        diag = (np.diag(np.diag(partial)))
        partial -= diag
        t_squared = partial*partial*(df / ((1.0 - partial) * (1.0 + partial)))
        prob = special.betainc(0.5*df, 0.5, df / (df + t_squared))
        np.fill_diagonal(partial, 1)
        np.fill_diagonal(prob, 0)
        return partial, prob
    else:
        return partial
     

def distance(matrix, metric="euclidean"):
    
    if metric=="euclidean":
        mat = spatial.distance.pdist(matrix.T)
    elif metric=="manhattan":
        mat = spatial.distance.pdist(matrix.T, metric='cityblock')
    
    return spatial.distance.squareform(mat)
    
    
def euclidean_distance(x,y):
    """Returns euclidean distance between two lists or numpy arrays"""
    x = np.array(x)
    y = np.array(y)
    return np.sqrt(sum((x-y)**2))


def manhattan_distance(x,y):
    """Returns manhattan distance between two lists or numpy arrays"""
    x = np.array(x)
    y = np.array(y)
    return sum(abs(x - y))


def calc_MI(X,Y,bins):
    """Returns Shannon's mutual information between two array-like objects"""
    c_XY = np.histogram2d(X,Y,bins)[0]
    c_X = np.histogram(X,bins)[0]
    c_Y = np.histogram(Y,bins)[0]

    H_X = shan_entropy(c_X)
    H_Y = shan_entropy(c_Y)
    H_XY = shan_entropy(c_XY)

    MI = H_X + H_Y - H_XY
    return MI

def shan_entropy(c):
    """Retuurns Shannon's information (entropy) for arrray-like object"""
    c_normalized = c / float(np.sum(c))
    c_normalized = c_normalized[np.nonzero(c_normalized)]
    H = -sum(c_normalized* np.log2(c_normalized))
    return H

def detrend(matrix)
    # remove global mean from the data:
    detrended = matrix - np.transpose(np.matlib.repmat(np.mean(matrix,axis=1), matrix.shape[1], 1))
    return detrended

def generate_fc(network_size, network_density, conn_weight):
    # create a synthetic matric of undirected functional connectivity (symmetric adjacency matrix):
    
    # convenient set of parameters:
    # network_size      = 10    number of nodes in the network
    # network_density   = 0.5   fraction of all possible connections
    # conn_weight       = 0.1   connectivity weight (between 0 and 1)


    # generate a vector for all the upper diagonal terms:
    num_comb = int(special.comb(network_size,2))
    
    upper_diag = np.zeros((num_comb,1))
    fc = np.zeros((network_size,network_size))
    
    fc_vec = np.random.rand(len(upper_diag),1)
    fc_vec[fc_vec > 1.0 - network_density] = 1
    fc_vec[fc_vec < 1.0 - network_density] = 0

    # fill in the upper triangular part of the FC matrix:
    inds_upperdiag = np.triu_indices(len(fc),1)
    fc[inds_upperdiag] = conn_weight*np.reshape(fc_vec, (num_comb,))
    inds_lowerdiag = np.triu_indices(len(fc),-1)

    # mirror upper triangular part to lower triangular:
    inds_lowerdiag = np.tril_indices(len(fc),-1)
    fc[inds_lowerdiag] = fc.T[inds_lowerdiag]

    # introduce self-inhibition on the diagonal:
    np.fill_diagonal(fc,-1)

    # indexes_conn = np.where(fc > 0)
    # indexes_noconn = np.where(fc == 0)
    
    return fc

    
def simulate_data(fcm, noise_level, T, dt):
    # fcm                       adjacency matrix generated with generate_fc function

    # convenient set of parameters:
    # noise_level       = 0.1   noise level in the system
    # T                 = 1000  length of the simulation [s]
    # dt                = 0.01  time step [s]
    
    network_size = fcm.shape[0]
    
    def deriv(x, t, adjacency):
        return np.dot(adjacency, x) + noise_inputs[int(t*100)-20,:] # + 0.0005*np.random.normal(0,1,(network_size,1))[:,0]

    adjacency = fcm

    time = np.linspace(0, T, int(np.floor(T/dt)) + 1)
    # initiate the variables at random:
    x0 = -0.5*np.ones((network_size,1))[:,0] + np.random.normal(0,1,(network_size,1))[:,0]
    # add noise to the system:
    noise_inputs = noise_level*(-0.5*np.ones((len(time),network_size)) + np.random.rand(len(time),network_size))
    synthetic_data = odeint(deriv, x0, time, args=(adjacency,))

    # shorten the data by the initial 25% samples:
    synthetic_data = synthetic_data[int(np.floor(len(time)*0.25)):,:]
    
    return synthetic_data