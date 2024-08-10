#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 22:01:20 2024

@author: sergiybegun
"""

# Clump Finding Problem with inclusion of overlapping patterns

def clump_finder(Genom_sequence, k_mer_length, L_sliding_window_length, t_required_occurence):
    """
    
    finding all the k-mer patterns that occured in all possible sliding windows
    of L_sliding_window_length size and their occurence is equal to t_required_occurence

    Parameters
    ----------
    Genom_sequence : str
        Genom sequence for the investigation.
    k_mer_length : int
        The length of pattern to search.
    L_sliding_window_length : int
        The length of sliding window within which to search.
    t_required_occurence : int
        The required number of occurence of k-mer inside current sliding window.

    Returns
    -------
    The list of all k-mer (k_mer_length) patterns with occurence equal to t_required_occurence
    within any possible sliding window of L_sliding_window_length size inside the Genom_sequence.
    For example,
    clump_finder("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4) = [CGACA, GAAGA]

    """
    
    clump_list = []
    
    length_of_Genom = len(Genom_sequence)
    
    if (length_of_Genom < k_mer_length) or (length_of_Genom < L_sliding_window_length):
        return []
    
    # cycle over all possible variants of sliding window
    for i in range(0,(length_of_Genom - L_sliding_window_length + 1)):
        
        dictionary_to_collect_preliminary_results = {}
        
        # cycle inside current sliding window
        for n in range(i,(i + L_sliding_window_length - k_mer_length)):
            current_pattern = Genom_sequence[n:(n + k_mer_length)]
            if current_pattern in dictionary_to_collect_preliminary_results.keys():
                dictionary_to_collect_preliminary_results[current_pattern] += 1
            else:
                dictionary_to_collect_preliminary_results[current_pattern] = 1
        
        # find if there are results within current window with occurence == t_required_occurence
        for curent_key in dictionary_to_collect_preliminary_results.keys():
            if (dictionary_to_collect_preliminary_results[curent_key] == t_required_occurence) and (curent_key not in clump_list):
                clump_list.append(curent_key)
    
    return clump_list


read_data_from_file = open("dataset_30274_5.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

string_with_numbers = read_strings_from_file[1].split()

clump_list_in_Genome = clump_finder(read_strings_from_file[0].strip(), int(string_with_numbers[0]), int(string_with_numbers[1]), int(string_with_numbers[2]))

output_file = open("output_sequence.txt", "w")

output_file.write(str(clump_list_in_Genome).replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

# real data example using the Genome of E. coli
# How many different 9-mers form (500,3)-clumps in the E. coli genome? (In other words, do not count a 9-mer more than once.)

read_data_from_file_1 = open("E_coli.txt", "r")

read_Genome_from_file_1 = read_data_from_file_1.read()

clump_list_in_Genome_1 = clump_finder(read_Genome_from_file_1.strip(), 9, 500, 3)

output_file_1 = open("output_E_coli_9_mers_500_3_clumps.txt", "w")

output_file_1.write(str(clump_list_in_Genome_1).replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file_1.close()

read_data_from_file_1.close()

number_of_different_E_coli_9_mers_500_3_clumps = len(clump_list_in_Genome_1)

print("number_of_different_E_coli_9_mers_500_3_clumps = ", number_of_different_E_coli_9_mers_500_3_clumps)
