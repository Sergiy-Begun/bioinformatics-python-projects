#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 08:14:24 2024

@author: sergiybegun
"""

def approximated_counts(maximum_mass):
    
    N_1 = 14712706211
    N_2 = 34544458837656
    
    m_1 = 1024
    m_2 = 1307
    
    C = ((float(N_2) / float(N_1)) ** (1.0 / (float(m_2) - float(m_1))))
    
    print("C = ", C)
    
    k = float(N_2) / (C ** float(m_2))
    
    print("k = ", k)
    
    Unknown_N_approximation = (k * (C ** maximum_mass))
    
    print("Unknown_N_approximation = ", Unknown_N_approximation)
    
    
approximated_counts(1473)