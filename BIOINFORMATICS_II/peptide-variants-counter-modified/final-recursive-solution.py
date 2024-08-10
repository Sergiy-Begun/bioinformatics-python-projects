#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:31:14 2024

@author: sergiybegun
"""

aminoacid_masses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]


def peptides(n, d):
    for m in aminoacid_masses:
        if n-m in d:
            print("n = ", n)
            print("m = ", m)
            print("n - m = ", (n - m))
            d[n] = d.get(n,0)+d[n-m]
            print("d[",n,"] = ", d[n])
    return d


def pep_counter(M):
    dicc = {0:1}
    mn = min(aminoacid_masses)
    print("mn = ", mn)
    for i in range(M-mn+1):
        j = i+mn
        peptides(j,dicc)
    return dicc

print(pep_counter(1290))