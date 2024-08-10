#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 11:09:47 2024

@author: sergiybegun
"""

import math

# calculate the entropy for each column in the genom matrix

def genom_matrix_entropy_by_columns(investigated_genom_matrix):
    """
    the estimation of entropy for columns of the genom matrix

    Parameters
    ----------
    investigated_genom_matrix : list of strings
        investigated genom matrix by strings of equal length.

    Returns
    -------
    the list of omputed values of entropy for each column of the string
    based on estimated probability values of occurrence 
    of nucleotides (A,T,C,and G) in the columns.

    """
    
    entropy = []
        
    genom_string_length = len(investigated_genom_matrix[0])
    
    number_of_strings = len(investigated_genom_matrix)
    
    for nucleotide_position in range(0,genom_string_length):
        
        current_counts = {
            "A": 0,
            "T": 0,
            "C": 0,
            "G": 0
            }
        
        current_entropy = 0.0
                
        for genom_string in investigated_genom_matrix:
            
            if (genom_string[nucleotide_position].lower() == "a"):
                current_counts["A"] += 1
            elif (genom_string[nucleotide_position].lower() == "t"):
                current_counts["T"] += 1
            elif (genom_string[nucleotide_position].lower() == "c"):
                current_counts["C"] += 1
            elif (genom_string[nucleotide_position].lower() == "g"):
                current_counts["G"] += 1
            else:
                raise ValueError("There is an error in genom string = " + str(genom_string))
        
        for nucleotide_key in current_counts.keys():
            current_entropy_parameter_for_nucleotide = 0.0
            current_probability_in_column_for_nucleotide = float(current_counts[nucleotide_key]) / float(number_of_strings)
            if (current_probability_in_column_for_nucleotide > 0):
                current_entropy_parameter_for_nucleotide = (-1.0) * (current_probability_in_column_for_nucleotide * math.log2(current_probability_in_column_for_nucleotide))
            current_entropy += current_entropy_parameter_for_nucleotide
        
        entropy.append(current_entropy)
    
    return entropy

read_data_from_file = open("input.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom_matrix = []
for i in range(0,len(read_strings_from_file)):
    investigated_genom_matrix.append(str(read_strings_from_file[i]).strip())

entropy_of_genom_matrix_strings = genom_matrix_entropy_by_columns(investigated_genom_matrix)

sum_of_entropies = 0.0
for entropy_per_column in entropy_of_genom_matrix_strings:
    sum_of_entropies += float(entropy_per_column)

output_file = open("output_entropy_by_column_list.txt", "w")

output_file.write(str(entropy_of_genom_matrix_strings).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.write("\nThe Sum of all entropy parameters for all columns = " + str(sum_of_entropies))

output_file.close()

read_data_from_file.close()
