#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:47:25 2024

@author: sergiybegun
"""

# finding most probable k-mer in a given genom string based on probbability matrix (4 x k)

def most_probable_k_mer_finder(investigated_genom_string: str, k_mer_length: int, probability_matrix):
    """
    

    Parameters
    ----------
    investigated_genom_string : str
        A genom string to find the most probable k-mer inside.
    k_mer_length : int
        A k-mer size.
    probability_matrix : TYPE
        An array of probabilities [4 x k].

    Returns
    -------
    A list of most probable k-mers in a given genom string.

    """
    
    most_probable_k_mer_list = []
    
    dictionary_to_convert_to_nucleotides_probability_matrix_strings = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3
        }
    
    dictionary_of_k_mers = {}
    
    length_of_genom_string = len(investigated_genom_string)
    
    if (k_mer_length > length_of_genom_string):
        raise ValueError("The k-mer length is greater than the genom length")
    
    for i in range(0,(length_of_genom_string - k_mer_length + 1)):
        
        current_moving_window_k_mer = investigated_genom_string[i:(i + k_mer_length)]
        
        current_probability = 1.0
        for j in range(0,k_mer_length):
            decoded_row_number_for_probability_matrix = dictionary_to_convert_to_nucleotides_probability_matrix_strings[current_moving_window_k_mer[j]]
            current_decoded_probability_for_nucleotide = probability_matrix[decoded_row_number_for_probability_matrix][j]
            
            current_probability *= current_decoded_probability_for_nucleotide
        
        if (current_moving_window_k_mer not in dictionary_of_k_mers.keys()):
            dictionary_of_k_mers[current_moving_window_k_mer] = current_probability
    
    max_probability = max(dictionary_of_k_mers.values())
    
    for current_key in dictionary_of_k_mers.keys():
        if (dictionary_of_k_mers[current_key] == max_probability):
            most_probable_k_mer_list.append(str(current_key))
    
    return most_probable_k_mer_list


read_data_from_file = open("dataset_30305_3.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom_input = str(read_strings_from_file[0]).strip()

k_mer_length_input =  int(str(read_strings_from_file[1]).strip())

probability_matrix_input = []

cur_pos = 0
for m in range(2,len(read_strings_from_file)):
    probability_matrix_input.append(read_strings_from_file[m].split())
    for l in range(0,k_mer_length_input):
        probability_matrix_input[cur_pos][l] = float(probability_matrix_input[cur_pos][l])
    cur_pos += 1

most_probable_k_mers = most_probable_k_mer_finder(investigated_genom_input,k_mer_length_input,probability_matrix_input)

output_file = open("most_probable_k_mer_list.txt", "w")

output_file.write(str(most_probable_k_mers).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()
