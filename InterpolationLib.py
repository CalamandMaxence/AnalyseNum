#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30

@author: Simone Deparis
"""

import numpy as np

# Defining the mxn Vandermonde matrix 
def VandermondeMatrix(x, m=-1):
    # Input
    # x : +1 array with interpolation nodes
    # m : degree of the polynomial. If negative or None, chooses m=size(x)-1 
    # Output
    # Matrix of Vandermonde of size m x n
                
    n = x.size - 1
    if (m is None or m < 0):
        m = n
                      
    # The following is a trick to put toghether the matrix A. Don't need to learn by heart ...
    # np.tile(x, (n+1, 1)) : repete n+1 times of the array x 
    #   .T : transpose
    # np.linspace(0,n,n+1) : [0,...,n+1]
    # np.tile : again: repete n+1 times the array  [0,...,n+1]
    # np.power : element by element power funcion

    A = np.power( np.tile(x, (m+1, 1)).T , np.tile(np.linspace(0,m,m+1), (n+1, 1)))

    return A




# Defining the Lagrange basis functions
def LagrangeBasis(t,k,z):
    # the input variables are:
    # t : the interpolatory points
    # k : which basis function
    # z : where to evaluate the function
    

    # careful, there are n+1 interpolation points!
    n = t.size - 1 

    # init result to one, of same type and size as z
    result = np.zeros_like(z) + 1

    # first few checks on k:
    if (type(k) != int) or (t.size < 1) or (k > n) or (k < 0):
        raise ValueError('Lagrange basis needs a positive integer k, smaller than the size of t')
        
    # loop on n to compute the product
    for j in range(0,n+1) :
        if (j == k) :
            continue
        if (t[k] == t[j]) :
            raise ValueError('Lagrange basis: all the interpolation points need to be distinct')

        result = result * (z - t[j]) / (t[k] - t[j])

    return result


# Lagrange Interpolation of data (t,y), evaluated at ordinate(s) z
def LagrangeInterpolation(t,y,z):
    '''
    the input variables are:
    t : the interpolatory points
    y : the corresponding data at the points t
    z : where to evaluate the function
    '''
    # {phi(t,k,.), k=0,...,n} is a basis of the polynomials of degree n
    # y represents the coordinates of the interpolating polynomial with respect to this basis.
    # Therefore LagrangePi(t,y,.) = y[0] phi(t,0,.) + ... + y[n] phi(t,n,.)
    
    # careful, there are n+1 basis functions!
    n = t.size - 1 

    # init result to zero, of same type and size as z
    result = np.zeros_like(z)
    
     # loop on n to compute the product
    for k in range(0,n+1) :
        result = result + y[k] * LagrangeBasis(t,k,z)

    return result  


# Piece-wise linear interpolation, equidistributed points
def PiecewiseLinearInterpolation(a,b,N,f,z):
    # the input variables are:
    # a,b : x[0] = a, x[n] = b
    # f : the corresponding data at the points x
    # z : where to evaluate the function

    def localLinearInterpolation (a, H, x, y, zk):
        # first find out in which interval we are. Easy here since equidistributed point
        i = int( (zk-a)//H )

        # if we are not in the interval, return constant function
        if x[i] == zk :
            return y[i]
        
        # use formula for local linear interpolation
        return y[i] + (y[i+1] - y[i])/(x[i+1] - x[i])* (zk - x[i]) 
    
    # careful, there are N+1 points 
    H = (b-a)/N
    if np.isclose(H,0) :
        raise ValueError('PiecewiseLinearInterpolation : a and b are too close')
    
    x = np.linspace(a,b,N+1)
    
    # if f is a function, we need to evaluate it at all the locations
    from types import FunctionType
    if isinstance(f, FunctionType) :
        y = f(x)
    else :
        y = f
        
    # if we want to evaluate the interpolating function in a single point
    if np.ndim(z) == 0 :
        return localLinearInterpolation (a, H, x, y, z)
    
    # if we are here, it means that we need to evaluate the interpolation on many points.
    # init result to zero, of same type and size as z
    result = np.zeros_like(z) 
     
    # The following part is not optimized
    for k in range(z.size) :
        result[k] = localLinearInterpolation (a, H, x, y, z[k])
   
    return result

# Piece-wise quadratic interpolation, equidistributed points
def PiecewiseQuadraticInterpolation(a,b,N,f,z):
    # the input variables are:
    # a,b : x[0] = a, x[n] = b
    # f : the corresponding data at the points x
    #     only works if f is a fonction
    # z : where to evaluate the function

    def localQuadraticInterpolation (a, H, x, f, zk):
        # first find out in which interval we are. Easy here since equidistributed point
        i = int( (zk-a)//H )

        # if we are not in the interval, return constant function
        if x[i] == zk :
            return f(x[i])
        
        t = np.array([ x[i], (x[i]+x[i+1])/2, x[i+1]] )
        
        # use formula for local quadratic interpolation
        return LagrangeInterpolation(t,f(t),zk)
        
    # careful, there are N+1 points 
    H = (b-a)/N
    if np.isclose(H,0) :
        raise ValueError('PiecewiseQuadraticInterpolation : a and b are too close')
    
    x = np.linspace(a,b,N+1)
    
    # if f is a function, we need to evaluate it at all the locations
    # from types import FunctionType
    # if ~isinstance(f, FunctionType) :
    #    raise ValueError('PiecewiseQuadraticInterpolation : works only with a given function')
        
    # if we want to evaluate the interpolating function in a single point
    if np.ndim(z) == 0 :
        return localQuadraticInterpolation (a, H, x, f, z)
    
    # if we are here, it means that we need to evaluate the interpolation on many points.
    # init result to zero, of same type and size as z
    result = np.zeros_like(z) 
     
    # The following part is not optimized
    for k in range(z.size) :

        result[k] = localQuadraticInterpolation (a, H, x, f, z[k])
   
    return result


   
    
    