#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 14:13:14 2024

@author: sergiybegun
"""

import time
import copy

t1 = time.time()

# find all contigs (genom fragments of maximum length with non-branching path based on de Bruijn graph)
# based on a given or generated from reads de Bruijn graph

# this task is required when there is no possibility to find unique Eulerian path

def contig_finder_function(de_bruijn_graph):
    """
    find all contigs based on de Bruijn graph (de_bruijn_graph)

    Parameters
    ----------
    de_bruijn_graph : dict
        The de Bruijn graph.

    Returns
    -------
    List of all contigs of maximum length.

    """
    
    list_of_all_contigs = []
    
    disbalance_viewer_dictionary = {}

    for cur_key in de_bruijn_graph.keys():
        
        current_de_bruijn_graph_string = de_bruijn_graph[cur_key]
        if (cur_key in disbalance_viewer_dictionary.keys()):
            disbalance_viewer_dictionary[cur_key][1].extend(current_de_bruijn_graph_string)
            for i in range (0,(len(current_de_bruijn_graph_string))):
                current_outlet_vertice = current_de_bruijn_graph_string[i]
                if (current_outlet_vertice in disbalance_viewer_dictionary.keys()):
                    disbalance_viewer_dictionary[current_outlet_vertice][0].append(cur_key)
                else:
                    disbalance_viewer_dictionary[current_outlet_vertice] = [[cur_key],[]]
        else:
            disbalance_viewer_dictionary[cur_key] = [[],current_de_bruijn_graph_string]
            for i in range (0,(len(current_de_bruijn_graph_string))):
                current_outlet_vertice = current_de_bruijn_graph_string[i]
                if (current_outlet_vertice in disbalance_viewer_dictionary.keys()):
                    disbalance_viewer_dictionary[current_outlet_vertice][0].append(cur_key)
                else:
                    disbalance_viewer_dictionary[current_outlet_vertice] = [[cur_key],[]]
    
    """
    for key_i in disbalance_viewer_dictionary.keys():
        inlet_counts = len(disbalance_viewer_dictionary[key_i][0])
        outlet_counts = len(disbalance_viewer_dictionary[key_i][1])
        in_out_difference = inlet_counts - outlet_counts
        if (in_out_difference < 0):
            inlet_disbalanced[key_i] = abs(in_out_difference)
        if (in_out_difference > 0):
            outlet_disbalanced[key_i] = in_out_difference
    """
    
    isolated_vertices_account = copy.deepcopy(list(de_bruijn_graph.keys()))
    
    list_of_de_gruijn_graph_vertices = copy.deepcopy(list(de_bruijn_graph.keys()))
    
    for i in range(0,len(list_of_de_gruijn_graph_vertices)):
        
        current_key_vertice = list_of_de_gruijn_graph_vertices[i]
        
        current_inlet_counts = len(disbalance_viewer_dictionary[current_key_vertice][0])
        current_outlet_counts = len(disbalance_viewer_dictionary[current_key_vertice][1])
        
        if (current_inlet_counts == 1) and (current_outlet_counts == 1):
            continue
        
        if (current_outlet_counts == 0):
            continue
        
        print("current starting vertice = ", current_key_vertice)
        
        if (current_key_vertice in isolated_vertices_account):
            isolated_vertices_account.remove(current_key_vertice)
        
        for current_next_vertice in disbalance_viewer_dictionary[current_key_vertice][1]:
            
            current_contig = str(current_key_vertice)
            
            next_vertice = current_next_vertice
            
            next_vertice_available = True
            
            while next_vertice_available == True:
                
                current_next_vertice_inlet_counts = len(disbalance_viewer_dictionary[next_vertice][0])
                current_next_vertice_outlet_counts = len(disbalance_viewer_dictionary[next_vertice][1])
                
                if (current_next_vertice_inlet_counts == 1) and (current_next_vertice_outlet_counts == 1):
                    
                    current_contig += " " + str(next_vertice)
                    
                    if (next_vertice in isolated_vertices_account):
                        isolated_vertices_account.remove(next_vertice)
                    
                    next_vertice = disbalance_viewer_dictionary[next_vertice][1][0]
                    
                    next_vertice_available = True
                else:
                    
                    current_contig += " " + str(next_vertice)
                    
                    if (next_vertice in isolated_vertices_account):
                        isolated_vertices_account.remove(next_vertice)
                    
                    list_of_all_contigs.append(current_contig)
                    
                    next_vertice_available = False
    
    print("len(isolated_vertices_account) = ", len(isolated_vertices_account))

    if (len(isolated_vertices_account) > 0):
        i = 0
        while (len(isolated_vertices_account) > 0):
            
            if (i >= len(isolated_vertices_account)):
                break
            
            current_next_vertice = isolated_vertices_account[i]
            
            next_vertice = disbalance_viewer_dictionary[current_next_vertice][1][0]
            
            current_contig = str(current_next_vertice) + " " + str(next_vertice)
            
            next_vertice_available = True
            
            k = 0
            
            while next_vertice_available == True:
                    
                if (next_vertice != current_next_vertice):
                    
                    if (k > 0):
                        current_contig += " " + str(next_vertice)
                    
                    if (next_vertice in isolated_vertices_account):
                        isolated_vertices_account.remove(next_vertice)
                    
                    next_vertice = disbalance_viewer_dictionary[next_vertice][1][0]
                    
                    next_vertice_available = True
                else:
                    
                    current_contig += " " + str(next_vertice)
                    
                    if (next_vertice in isolated_vertices_account):
                        isolated_vertices_account.remove(next_vertice)
                    
                    if (current_next_vertice in isolated_vertices_account):
                        isolated_vertices_account.remove(current_next_vertice)
                    
                    list_of_all_contigs.append(current_contig)
                    
                    next_vertice_available = False
                
                k += 1
            
            i += 1
            
    
    return sorted(list_of_all_contigs)

read_data_from_file = open("dataset_30209_2.txt", "r")    

read_strings_from_file = read_data_from_file.readlines()

vertices_dictionary = {}

for graph_string in read_strings_from_file:
    current_input_string = graph_string.split()
    current_first_element = str(current_input_string[0]).replace(":", "")
    vertices_dictionary[current_first_element] = []
    for m in range(1,len(current_input_string)):
        vertices_dictionary[current_first_element].append(str(current_input_string[m]).strip())

print("vertices_dictionary = ", vertices_dictionary)

contig_list_result = contig_finder_function(vertices_dictionary)

print("contig_list_result = ", contig_list_result)

output_file = open("output_contig_list.txt", "w")

for list_element in contig_list_result:
    output_file.write(str(list_element).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))
    output_file.write("\n")

output_file.close()

read_data_from_file.close()

print("running time = ", (time.time() - t1))
