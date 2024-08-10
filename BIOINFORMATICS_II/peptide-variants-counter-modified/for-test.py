#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 17:14:40 2024

@author: sergiybegun
"""

def for_test(n):    
    for i_1 in range(0,10):
        if n < 1:
            break
        if (n >= 1):
            print("n = ", n, " i_1 = ", i_1)
        for i_2 in range(0,10):
            if n < 2:
                break
            if (n >= 2):
                print("n = ", n, " i_1(",i_1, " + i_2(",i_2, " = ", (i_1 + i_2))
            for i_3 in range(0,10):
                if n < 3:
                    break
                if (n >= 3):
                    print("n = ", n, " i_1(",i_1, ") + i_2(",i_2, ") + i_3(",i_3,") = ", (i_1 + i_2 + i_3))
            for i_4 in range(0,10):
                if n < 4:
                    break
                if (n >= 4):
                    print("n = ", n, " i_1(",i_1, ") + i_2(",i_2, ") + i_3(",i_3, ") + i_4(",i_4,") = ", (i_1 + i_2 + i_3 + i_4))

for_test(2)