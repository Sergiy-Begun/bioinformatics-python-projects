#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 16:00:52 2024

@author: sergiybegun
"""

# finding all positions (starting points) in the genom sequence 
# where HammingDistance(Pattern, Pattern') ≤ d

def Hamming_distance_calculation(genom_sequence_1, genom_sequence_2):
    """
    Calculating the Hamming distance between two genome sequence of equal length
    as a count of nucleotide mistmatches at the equal positions along the strandds.

    Parameters
    ----------
    genom_sequence_1 : str
        First genome fragment.
    genom_sequence_2 : str
        Second genome fragment.

    Returns
    -------
    Number of nucleotide mistmatches in terms of Hamming distance.
    For example, Hamming_distance_calculation("GGGCCGTTGGT", "GGACCGTTGAC") = 3

    """
    
    Hamming_distance_measure = 0
    
    length_gen_seq_1 = len(genom_sequence_1)
    length_gen_seq_2 = len(genom_sequence_2)
    
    if (length_gen_seq_1 != length_gen_seq_2):
        return -1
    
    for i in range(0,length_gen_seq_1):
        if (genom_sequence_1[i] != genom_sequence_2[i]):
            Hamming_distance_measure += 1
    
    return Hamming_distance_measure

def find_possible_dnaa_candidates(pattern_to_find, input_genome_sequence, allowed_Hamming_distance):
    """
    # finding all positions (starting points) in the input_genome_sequence 
    # where HammingDistance(pattern_to_find, Pattern') ≤ allowed_Hamming_distance

    Parameters
    ----------
    pattern_to_find : str
        An approximate pattern-candidate to find.
    input_genome_sequence : str
        A fragment of genome strand to investigate.
    allowed_Hamming_distance : int
        Allowed value of Hamming distance between pattern_to_find and fragments 
        investigated within sliding window of length equal to len(pattern_to_find).

    Returns
    -------
    The list of all starting positions off aproximate-candidates patterns.
    For example, find_possible_dnaa_candidates("ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3) = [6, 7, 26, 27]

    """
    
    list_to_return = []
    
    length_of_genom = len(input_genome_sequence)
    length_of_pattern = len(pattern_to_find)
    
    if (length_of_genom < length_of_pattern) or (length_of_genom < allowed_Hamming_distance):
        return []
    
    for i in range(0,(length_of_genom - length_of_pattern + 1)):
        current_pattern = input_genome_sequence[i:(i + length_of_pattern)]
        if (Hamming_distance_calculation(pattern_to_find, current_pattern) <= allowed_Hamming_distance):
            list_to_return.append(i)
    
    return list_to_return

read_data_from_file = open("dataset_30278_6.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_pattern = str(read_strings_from_file[0]).strip()

investigated_genom = str(read_strings_from_file[1]).strip()

allowed_distance = int(str(read_strings_from_file[2]).strip())

list_of_starting_positions = find_possible_dnaa_candidates(investigated_pattern, investigated_genom, allowed_distance)

output_file = open("output_dnaa_starting_positions_list.txt", "w")

output_file.write(str(list_of_starting_positions).replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

# count the number of approximate dnaa candidates

number_of_dnaa_candidates = len(find_possible_dnaa_candidates("AAAAA", "AACAAGCTGATAAACATTTAAAGAG", 2))

print("number_of_dnaa_candidates = ", number_of_dnaa_candidates)

number_of_dnaa_candidates_from_main_output = len(list_of_starting_positions)

print("number_of_dnaa_candidates_from_main_output  = ", number_of_dnaa_candidates_from_main_output )

