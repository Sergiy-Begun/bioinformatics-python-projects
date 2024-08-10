#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 20:32:42 2024

@author: sergiybegun
"""

# reconstruct genom strand step-by-step from the beginning to the end from 
# the pieces of k-mers of equivalent length given in a correct numbering

def simple_genom_sewing(list_of_k_mers):
    """
    
    sewing of k-mers into a genom strand

    Parameters
    ----------
    list_of_k_mers : list of str
        The list of k-mers to sew into the genom strand.

    Returns
    -------
    Sewed genom strand.
    For example, simple_genom_sewing(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"]) = ACCGAAGCT

    """
    
    sewed_genom_strand = ""
    
    sewed_genom_strand += list_of_k_mers[0]
    
    for i in range(1,len(list_of_k_mers)):
        sewed_genom_strand += list_of_k_mers[i][-1]
    
    return sewed_genom_strand

read_data_from_file = open("dataset_30182_3.txt", "r")

read_strings_from_file = read_data_from_file.read().split()

investigated_genom = []

for m in range(0,len(read_strings_from_file)):
    investigated_genom.append(str(read_strings_from_file[m]).strip())

genom_strand = simple_genom_sewing(investigated_genom)

output_file = open("output_genom_strand.txt", "w")

output_file.write(str(genom_strand).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()
