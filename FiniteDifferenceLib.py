#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30

@author: Simone Deparis
"""

import numpy as np

# This function just provides a line to compare the convergence in a log-log plot.
def sampleConvergence(h,order,ref) :
    # computes a pseudo order of convergence on h
    # ref is a reference maximum value 
    return h**order * (ref/np.max(np.abs(h))**order)


   
    
    
