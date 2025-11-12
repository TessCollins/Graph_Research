import numpy as np
from sage.modules.misc import gram_schmidt
from sage.all import *

##################################################

def lib_make_uniform_rand_vec(dimension):
    """
    Makes a uniformly random vector.
    
    Args:
    - Desired dimension
    
    Returns:
    - Uniformly random vector of that dimension 
    """
    xs = np.random.normal(0,1,dimension) #normally random x's
    norm = sum(xs**2)**0.5 #and then normalize
    return vector(xs/norm)

##################################################

def lib_make_flat_vec(dim):
    """
    Makes a flat vector.
    
    Args:
    - Desired dimension
    
    Returns:
    - Vector of that dimension with entries of modulus 1
    """
    theta = np.random.uniform(0,2*pi,dim)    
    vec = []
    for j in range(dim):
        vec.append(math.sin(theta[j])+math.cos(theta[j])*I)

    return vec

##################################################

def lib_make_Fourier_Matrix(d):
    """
    Makes the Fourier matrix associated with a specified dimension

    Args:
    - Dimension

    Returns:
    - Fourier matrix of that dimension
    """

    rs=[]
    for k in range(d):
        rw = [None for i in range(d)]
        for j in range(d):
            rw[j] = exp(2*math.pi*I*k*j/d)
        rs.append(rw)
    F = rs
    for k in range(d):
        for j in range(d):
            if abs(F[k][j].real()) < 1e-5:
                F[k][j]=F[k][j].imag()*I
            if abs(F[k][j].imag()) < 1e-5:
                F[k][j]=F[k][j].real()
    return (1/sqrt(d))*Matrix(F)

##################################################

def lib_make_standard_basis_vectors(d):
    """
    Will give a list of the standard basis vectors in dimension d.
    
    Args:
    - Dimension d
    
    Returns:
    - List of standard basis vectors, each vector as a list
    """
    veclist = []
    for i in range(d):
        vec = []
        for j in range(d):
            if i == j:
                vec.append(1)
            else:
                vec.append(0)
        veclist.append(Matrix(vec).transpose())
    return veclist    

##################################################