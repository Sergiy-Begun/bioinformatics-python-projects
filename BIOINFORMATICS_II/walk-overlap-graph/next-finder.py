#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 21:52:09 2024

@author: sergiybegun
"""

# give for every k-mer in the list the list of k-mers, where suffix(i)=prefix(j)

def next_finder(list_of_k_mers):
    """
    A dictionary with list of all k-mers, where suffix(i)=prefix(j) for a given k-mer.

    Parameters
    ----------
    list_of_k_mers : list of str
        A list of investigated k-mers.

    Returns
    -------
     A dictionary with list of all k-mers, where suffix(i)=prefix(j) for a given k-mer.

    """
    
    vertices_dictionary = {}
    
    # formation of tuples (prefix,suffix)
    
    for k_mer in list_of_k_mers:
        if (k_mer not in vertices_dictionary.keys()):
            vertices_dictionary[k_mer] = [(k_mer[0:-1],k_mer[1:len(k_mer)]),[]]
    
    for key_i in vertices_dictionary.keys():
        for key_j in vertices_dictionary.keys():
            if (key_i == key_j):
                continue
            
            cur_suf_i = vertices_dictionary[key_i][0][1]
            cur_pref_j = vertices_dictionary[key_j][0][0]
            
            if (cur_suf_i == cur_pref_j):
                vertices_dictionary[key_i][1].append(key_j)
    
    return vertices_dictionary

read_data_from_file = open("dataset_30182_10.txt", "r")

read_strings_from_file = read_data_from_file.read().split()

investigated_genom = []

for m in range(0,len(read_strings_from_file)):
    investigated_genom.append(str(read_strings_from_file[m]).strip())

vertices_dictionary_output = next_finder(investigated_genom)

output_file = open("output_vertices_dictionary.txt", "w")

for cur_vertice in vertices_dictionary_output.keys():
    if (len(vertices_dictionary_output[cur_vertice][1]) == 0):
        continue
    output_file.write((str(cur_vertice) + ": "))
    output_file.write((str(vertices_dictionary_output[cur_vertice][1]).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", "") + "\n"))

output_file.close()

read_data_from_file.close()
