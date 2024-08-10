#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 22:05:31 2024

@author: sergiybegun
"""

# build De Bruijn Graph with vertices equal to prefixes, directions equal to suffixes,
# and edges equal to k-mers for a given length of k-mers


def graph_vertices_based_on_prefixes(list_of_k_mers):
    """
    returns the dictionary of graph of the k-mers

    Parameters
    ----------
    list_of_k_mers : list of str
        investigated genom k-mers of equal length.

    Returns
    -------
    A dictionary of graph vertices with lists of the connections and edges
    {prefix: [(edge,suffix)1,(edge,suffix)2...,(edge,suffix)n]}

    """
    
    vertices_dictionary = {}
    
    # formation of tuples (edge,suffix)
    
    for k_mer in list_of_k_mers:
        k_mer_prefix = k_mer[0:-1]
        k_mer_suffix = k_mer[1:len(k_mer)]
        if (k_mer_prefix in vertices_dictionary.keys()):
            vertices_dictionary[k_mer_prefix].append((k_mer,k_mer_suffix))
            
        if (k_mer_prefix not in vertices_dictionary.keys()):
            vertices_dictionary[k_mer_prefix] = [(k_mer,k_mer_suffix)]
        
    return vertices_dictionary


read_data_from_file = open("dataset_30184_8.txt", "r")

read_strings_from_file = read_data_from_file.read()

k_mers_list = str(read_strings_from_file).split()

for m in range(0,len(k_mers_list)):
    k_mers_list[m] = str(k_mers_list[m]).strip()

vertices_dictionary_output = graph_vertices_based_on_prefixes(k_mers_list)

output_file = open("output_prefix_based_vertices_dictionary.txt", "w")

for cur_vertice in vertices_dictionary_output.keys():
    output_file.write((str(cur_vertice) + ": "))
    for m in range(0,len(vertices_dictionary_output[cur_vertice])):
        output_file.write(str(vertices_dictionary_output[cur_vertice][m][1]).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))
        if (m != (len(vertices_dictionary_output[cur_vertice]) - 1)):
            output_file.write(" ")
    output_file.write("\n")

output_file.close()

read_data_from_file.close()
