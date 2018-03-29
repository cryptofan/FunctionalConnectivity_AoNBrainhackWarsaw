# This file is intended to work as a library.
# you can
# import fc_methods as fc
# and then e.g.:
# fc.manhattan_distance(x,y)
# etc.
# Happy coding! :)


# Import dependencies
from scipy import special, spatial#, stats, linalg
import numpy as np

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
