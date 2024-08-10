#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 17:02:24 2024

@author: sergiybegun
"""

import copy
import random
import time

t1 = time.time()


def de_broijn_graph_balancer(potentially_unbalanced_graph):
    disbalance_viewer_dictionary = {}
    inlet_disbalanced = {}
    outlet_disbalanced = {}
    for cur_key in potentially_unbalanced_graph.keys():
        
        current_de_bruijn_graph_string = potentially_unbalanced_graph[cur_key]
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
        
    for key_i in disbalance_viewer_dictionary.keys():
        inlet_counts = len(disbalance_viewer_dictionary[key_i][0])
        outlet_counts = len(disbalance_viewer_dictionary[key_i][1])
        in_out_difference = inlet_counts - outlet_counts
        if (in_out_difference < 0):
            inlet_disbalanced[key_i] = abs(in_out_difference)
        if (in_out_difference > 0):
            outlet_disbalanced[key_i] = in_out_difference
    
    print("inlet_disbalanced = ", inlet_disbalanced)
    print("outlet_disbalanced = ", outlet_disbalanced) 
    print("disbalance_viewer_dictionary = ", disbalance_viewer_dictionary)
    if ((max(inlet_disbalanced.values()) > 1) or (max(outlet_disbalanced.values()) > 1) or (len(inlet_disbalanced.keys()) > 1) or (len(outlet_disbalanced.keys()) > 1)):
        raise ValueError("Algorithm shoud be enhanced for more than one disbalanced pair")
    
    artificial_vertice = list(outlet_disbalanced.keys())[0]
    artificial_outlet = list(inlet_disbalanced.keys())[0]
    if (artificial_vertice in potentially_unbalanced_graph.keys()):
        potentially_unbalanced_graph[artificial_vertice].append(artificial_outlet)
    else:
        potentially_unbalanced_graph[artificial_vertice] = [artificial_outlet]
    
    edge_for_deletion_in_future_eulerian_cycle_to_form_eulerian_path = str(list(outlet_disbalanced.keys())[0]) + str(list(inlet_disbalanced.keys())[0])
    
    return (potentially_unbalanced_graph, edge_for_deletion_in_future_eulerian_cycle_to_form_eulerian_path)
        

def eulerian_cycle_finder(vertices_dictionary_input):
    full_list_of_edges = []
    starting_positions_with_edges = {}
    
    for cur_key in vertices_dictionary_input.keys():
        current_de_bruijn_graph_string = vertices_dictionary_input[cur_key]
        starting_positions_with_edges[cur_key] = []
        for i in range(0,(len(current_de_bruijn_graph_string))):
            new_edge = str(cur_key) + " " + current_de_bruijn_graph_string[i]
            starting_positions_with_edges[cur_key].append(new_edge)
            #if new_edge not in full_list_of_edges:
            full_list_of_edges.append(new_edge)
    
    unused_edges_full_list = copy.deepcopy(full_list_of_edges)
    current_cycle = {}
    
    list_of_starting_keys = list(starting_positions_with_edges.keys())
    random_position_in_starting_list = random.randint(0,(len(list_of_starting_keys) - 1))
    
    starting_key = list_of_starting_keys[random_position_in_starting_list]
    
    first_list_of_edges = starting_positions_with_edges[starting_key]
    random_edge_number = random.randint(0, (len(first_list_of_edges) - 1))
    first_vertice_unused_edges = copy.deepcopy(starting_positions_with_edges[starting_key])
    first_vertice_unused_edges.remove(first_list_of_edges[random_edge_number])
    next_element_id_of_the_cycle = str(first_list_of_edges[random_edge_number]).split()[1]
    current_cycle[starting_key] = [starting_positions_with_edges[starting_key],first_vertice_unused_edges]
    unused_edges_full_list.remove(first_list_of_edges[random_edge_number])
    current_version_of_full_path = [starting_key,next_element_id_of_the_cycle]
    
    next_element_of_cycle_could_be_generated = True
        
    while next_element_of_cycle_could_be_generated == True:
        #print("len(full_list_of_edges) = ", len(full_list_of_edges))
        #print("len(current_version_of_full_path) = ", len(current_version_of_full_path))

        if (len(unused_edges_full_list) == 0):
            return current_version_of_full_path

        if (next_element_id_of_the_cycle in current_cycle.keys()):
            next_element_id_of_the_cycle_before = next_element_id_of_the_cycle
            existing_current_cycle_unused_edges = current_cycle[next_element_id_of_the_cycle_before][1]
            length_of_existing_current_cycle_unused_edges = len(existing_current_cycle_unused_edges)
            if (length_of_existing_current_cycle_unused_edges > 0):
                random_existing_edge_number = random.randint(0, (len(existing_current_cycle_unused_edges) - 1))
                randomly_chosen_edge_for_existing_vertice = existing_current_cycle_unused_edges[random_existing_edge_number]
                existing_current_cycle_unused_edges.remove(randomly_chosen_edge_for_existing_vertice)
                next_element_id_of_the_cycle_after = str(randomly_chosen_edge_for_existing_vertice).split()[1]
                current_cycle[next_element_id_of_the_cycle_before][1] = copy.deepcopy(existing_current_cycle_unused_edges)
                unused_edges_full_list.remove(randomly_chosen_edge_for_existing_vertice)
                next_element_id_of_the_cycle = next_element_id_of_the_cycle_after
                current_version_of_full_path.append(next_element_id_of_the_cycle_after)
                next_element_of_cycle_could_be_generated = True
            else:
                next_element_of_cycle_could_be_generated = False
                copy_current_version_of_full_path = copy.deepcopy(current_version_of_full_path)
                current_version_of_full_path = []
                
                unused_edge_position = 0
                
                for cur_i in range(0,len(copy_current_version_of_full_path)):
                    #print("cur_i = ", cur_i)
                    existing_path_key_element = copy_current_version_of_full_path[cur_i]
                    #print("existing_path_key_element = ", existing_path_key_element)
                    according_current_cycle_vertice_unused_edges = current_cycle[existing_path_key_element][1]
                    #print("copy_current_version_of_full_path = ", str(copy_current_version_of_full_path))
                    #print("according_current_cycle_vertice_unused_edges = ", str(according_current_cycle_vertice_unused_edges))
                    if (len(according_current_cycle_vertice_unused_edges) > 0):
                        new_starting_element_of_the_path = existing_path_key_element
                        unused_edge_position = cur_i
                        current_version_of_full_path.extend(copy_current_version_of_full_path[unused_edge_position:(len(copy_current_version_of_full_path) + 1)])
                        current_version_of_full_path.extend(copy_current_version_of_full_path[1:unused_edge_position])
                        current_version_of_full_path.append(new_starting_element_of_the_path)
                        next_element_of_cycle_could_be_generated = True
                        next_element_id_of_the_cycle = new_starting_element_of_the_path
                        break
        else:
            next_element_id_of_the_cycle_before = next_element_id_of_the_cycle
            new_list_of_edges = starting_positions_with_edges[next_element_id_of_the_cycle_before]
            random_new_edge_number = random.randint(0, (len(new_list_of_edges) - 1))
            new_vertice_unused_edges = copy.deepcopy(starting_positions_with_edges[next_element_id_of_the_cycle_before])
            randomly_chosen_edge = new_list_of_edges[random_new_edge_number]
            new_vertice_unused_edges.remove(randomly_chosen_edge)
            next_element_id_of_the_cycle_after = str(randomly_chosen_edge).split()[1]
            current_cycle[next_element_id_of_the_cycle_before] = [new_list_of_edges,new_vertice_unused_edges]
            unused_edges_full_list.remove(randomly_chosen_edge)
            next_element_id_of_the_cycle = next_element_id_of_the_cycle_after
            current_version_of_full_path.append(next_element_id_of_the_cycle_after)
            next_element_of_cycle_could_be_generated = True
    return current_version_of_full_path

def eulerian_path_formation_from_eulerian_cycle(eulerian_cycle_input, artificially_injected_edge):
    eulerian_path_result = []
    length_of_eulerian_cycle = len(eulerian_cycle_input)
    for i in range(0,(length_of_eulerian_cycle - 1)):
        searching_for_artificially_injected_edge = str(eulerian_cycle_input[i]) + str(eulerian_cycle_input[i + 1])
        if (searching_for_artificially_injected_edge == artificially_injected_edge):
            eulerian_path_result.extend(eulerian_cycle_input[(i + 1):(length_of_eulerian_cycle + 1)])
            eulerian_path_result.extend(eulerian_cycle_input[1:(i + 1)])
    
    return eulerian_path_result

read_data_from_file = open("test_input.txt", "r")    

read_strings_from_file = read_data_from_file.readlines()

vertices_dictionary = {}

for graph_string in read_strings_from_file:
    current_input_string = graph_string.split()
    current_first_element = str(current_input_string[0]).replace(":", "")
    vertices_dictionary[current_first_element] = []
    for m in range(1,len(current_input_string)):
        vertices_dictionary[current_first_element].append(str(current_input_string[m]).strip())

print("original graph = ", vertices_dictionary)
balaced_graph_result = de_broijn_graph_balancer(vertices_dictionary)

balanced_graph = balaced_graph_result[0]
print("balanced_graph = ", balanced_graph)

edge_for_deletion_which_contains_begin_and_end = balaced_graph_result[1]

eulerian_cycle = eulerian_cycle_finder(balanced_graph)

#print("eulerian_cycle = ", eulerian_cycle)

#print("edge_for_deletion_which_contains_begin_and_end = ", edge_for_deletion_which_contains_begin_and_end)

eulerian_path_output = eulerian_path_formation_from_eulerian_cycle(eulerian_cycle,edge_for_deletion_which_contains_begin_and_end)

#print("eulerian_path_output = ", eulerian_path_output)

output_file = open("output_eulerian_path.txt", "w")

output_file.write(str(eulerian_path_output).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

print("running time = ", (time.time() - t1))
