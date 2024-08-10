#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 20:21:39 2024

@author: sergiybegun
"""

# finding all the starting positions in Genome where Pattern appears as a substring

def starting_positions_finder(input_Genome, pattern_to_find):
    """
    finding all the starting positions in Genome where Pattern appears as a substring

    Parameters
    ----------
    input_Genome : str
        input sequence of Genome.
    
    pattern_to_find : str
        pattern sequence to find in Genome.

    Returns
    -------
    An array (list) of all starting positions of pattern_to_find in input_Genome.
    For example, starting_positions_finder("GATATATGCATATACTT", "ATAT") = [1, 3, 9]

    """
    
    list_of_starting_points = []
    
    length_of_the_Genome = len(input_Genome)
    
    length_of_the_pattern = len(pattern_to_find)
    
    if (length_of_the_pattern > length_of_the_Genome):
        return []
    
    for i in range(0,(length_of_the_Genome - length_of_the_pattern + 1)):
        if input_Genome[i:(i + length_of_the_pattern)].lower() == pattern_to_find.lower():
            list_of_starting_points.append(i)
    
    return list_of_starting_points

read_data_from_file = open("dataset_30273_5.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

starting_positions_list = starting_positions_finder(read_strings_from_file[1].strip(), read_strings_from_file[0].strip())

output_file = open("output_sequence", "w")

output_file.write(str(starting_positions_list).replace(", ", " ").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

# real data example using the Genome of Vibrio_cholerae

read_data_from_file_1 = open("Vibrio_cholerae.txt", "r")

read_Genome_from_file_1 = read_data_from_file_1.read()

pattern_to_find_real = "CTTGATCAT"

starting_positions_list_1 = starting_positions_finder(read_Genome_from_file_1.strip(), pattern_to_find_real)

output_file_1 = open("output_sequence_Vibrio_cholerae.txt", "w")

output_file_1.write(str(starting_positions_list_1).replace(", ", " ").replace("[", "").replace("]", ""))

output_file_1.close()

read_data_from_file_1.close()