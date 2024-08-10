#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 14:21:52 2024

@author: sergiybegun
"""

# build De Bruijn Graph with vertices equal to prefixes, directions equal to suffixes,
# and edges equal to k-mers for a given length of k-mers


def all_possible_k_mers(genom_strand: str, k_mer_length: int):
    """
    
    give a list of all possible k-mer fragments for a given genom strand

    Parameters
    ----------
    genom_strand : str
        investigated genom strand.
    k_mer_length : int
        k-mer length.

    Returns
    -------
    a list of all possible k-mer fragments for a given genom strand.
    For example, all_possible_k_mers("CAATCCAAC",5) = ["CAATC", "AATCC", "ATCCA", "TCCAA", "CCAAC"]

    """
    
    list_of_k_mers = []
    
    genom_strand_length = len(genom_strand)
    
    if (genom_strand_length < k_mer_length):
        raise ValueError("Genom strand length is shorter than the k-mer length")
    
    for i in range(0,(genom_strand_length - k_mer_length + 1)):
        list_of_k_mers.append(genom_strand[i:(i + k_mer_length)])
    
    return list_of_k_mers

def graph_vertices_based_on_prefixes(genom_strand: str, k_mer_length: int):
    """
    returns the dictionary of graph of the k-mers

    Parameters
    ----------
    genom_strand : str
        investigated genom strand.
    k_mer_length : int
        k-mer length.

    Returns
    -------
    A dictionary of graph vertices with lists of the connections and edges
    {prefix: [(edge,suffix)1,(edge,suffix)2...,(edge,suffix)n]}

    """
    
    vertices_dictionary = {}
    
    length_of_the_genom_strand = len(genom_strand)
    
    if (length_of_the_genom_strand < k_mer_length):
        raise ValueError("Genom strand's length is shorter than the k-mer length")
    
    all_possible_k_mers_list = all_possible_k_mers(genom_strand,k_mer_length)
    
    # formation of tuples (edge,suffix)
    
    for k_mer in all_possible_k_mers_list:
        k_mer_prefix = k_mer[0:-1]
        k_mer_suffix = k_mer[1:len(k_mer)]
        if (k_mer_prefix in vertices_dictionary.keys()):
            vertices_dictionary[k_mer_prefix].append((k_mer,k_mer_suffix))
            
        if (k_mer_prefix not in vertices_dictionary.keys()):
            vertices_dictionary[k_mer_prefix] = [(k_mer,k_mer_suffix)]
        
    return vertices_dictionary


read_data_from_file = open("dataset_30183_6.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom = str(read_strings_from_file[1]).strip()

k_mer_length_input = int(str(read_strings_from_file[0]).strip())

vertices_dictionary_output = graph_vertices_based_on_prefixes(investigated_genom,k_mer_length_input)

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
