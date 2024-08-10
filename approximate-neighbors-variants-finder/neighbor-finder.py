#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 23:15:59 2024

@author: sergiybegun
"""

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
        return ("Error = length mismatch")
    
    for i in range(0,length_gen_seq_1):
        if (genom_sequence_1[i] != genom_sequence_2[i]):
            Hamming_distance_measure += 1
    
    return Hamming_distance_measure


def reverse_complement(input_dna_sequence):
    """
    creating the reverse complement DNA sequence

    Parameters
    ----------
    input_dna_sequence : str
        input DNA sequence.

    Returns
    -------
    The reverse complement DNA sequence.
    Example reverse_complement("AAAACCCGGT") = "ACCGGGTTTT"

    """
    
    # our decoding dictionary
    decoding_dictionary = {
        "A": "T",
        "a": "T",
        "C": "G",
        "c": "G",
        "T": "A",
        "t": "A",
        "G": "C",
        "g": "C"
        }
    
    output_dna_sequence = ""
    
    # creating the direct complementary sequence
    for i in range(0,len(input_dna_sequence)):
        current_symbol = input_dna_sequence[i]
        
        if current_symbol not in decoding_dictionary.keys():
            return ("Error = " + current_symbol + " i = ", i)
        
        output_dna_sequence += decoding_dictionary[current_symbol]
        
    # creating reverse sequense DNA
    
    output_dna_sequence = output_dna_sequence[::-1]
    
    return output_dna_sequence


def create_neighbors(investigated_pattern, allowed_hamming_distance,strict_equality):
    """
    Generating from investigated_pattern all possible neighbors
    in terms of allowed Hamming distance <= (or ==) allowed_hamming_distance
    as two lists of patterns (first list), and a complementaries to that patterns (second list).

    Parameters
    ----------
    investigated_pattern : str
        input k-mer pattern for investigation.
    allowed_hamming_distance : int
        Allowed Hamming distance for such approximate k-mers.
    strict_equality : boolean
        "<=" if strict_equality is False, and "==" if strict_equality is True.

    Returns
    -------
    The list of all possible neighbors
    in terms of allowed Hamming distance <= (or ==) allowed_hamming_distance
    as two lists of patterns (first list), and a complementaries to that patterns (second list).
    The investigated_pattern is the first item in the first list.
    The complementary to the investigated_pattern is the first item in the second list.

    """
    
    list_of_nucleotide_variants = ["A", "T", "C", "G"]
    
    # The input_pattern is the first item in the first list.
    dictionary_of_neighbors = {investigated_pattern: 0}
    # The complementaru to input_pattern is the first item in the second list.
    dictionary_of_neighbors_complementary = {reverse_complement(investigated_pattern): 0}
    
    length_of_input_pattern = len(investigated_pattern)
    
    if (length_of_input_pattern < allowed_hamming_distance):
        return ([], [])
        
    if (allowed_hamming_distance == 0):
        return (dictionary_of_neighbors.keys(),dictionary_of_neighbors_complementary.keys())
    
    current_hamming_distance_level = 0
    
    while current_hamming_distance_level < allowed_hamming_distance:
        
        current_hamming_distance_level += 1
        
        current_list_of_direct_dictionary_keys = []
        for i_key in dictionary_of_neighbors.keys():
            current_list_of_direct_dictionary_keys.append(i_key)
                
        for current_direct_variant in current_list_of_direct_dictionary_keys:
            if (dictionary_of_neighbors[current_direct_variant] != (current_hamming_distance_level - 1)):
                continue

            for j in range(0, length_of_input_pattern):
                for nucleotide in list_of_nucleotide_variants:
                    current_direct_nucleotide = str(current_direct_variant)[j]
                    
                    if (current_direct_nucleotide != nucleotide) and (current_direct_nucleotide == investigated_pattern[j]):
                        new_direct_neighbor = (current_direct_variant[0:j] + nucleotide + current_direct_variant[(j + 1):length_of_input_pattern])
                        
                        if new_direct_neighbor not in dictionary_of_neighbors.keys():
                            dictionary_of_neighbors[new_direct_neighbor] = current_hamming_distance_level
                            new_complementary_neighbor = reverse_complement(new_direct_neighbor)
                            dictionary_of_neighbors_complementary[new_complementary_neighbor] = current_hamming_distance_level
    
    list_of_direct_pattens = []
    list_of_complementary_patterns = []
    
    if strict_equality == False:
        for direct_key in dictionary_of_neighbors.keys():
            list_of_direct_pattens.append(str(direct_key))
                
        for complementary_key in dictionary_of_neighbors_complementary.keys():
            list_of_complementary_patterns.append(str(complementary_key))
    
    if strict_equality == True:

        for direct_key in dictionary_of_neighbors.keys():
            if (dictionary_of_neighbors[direct_key] == allowed_hamming_distance) or (dictionary_of_neighbors[direct_key] == 0):
                list_of_direct_pattens.append(str(direct_key))
                
        for complementary_key in dictionary_of_neighbors_complementary.keys():
            if (dictionary_of_neighbors_complementary[complementary_key] == allowed_hamming_distance) or (dictionary_of_neighbors_complementary[complementary_key] == 0):
                list_of_complementary_patterns.append(str(complementary_key))
                
    tuple_to_return = (list_of_direct_pattens,list_of_complementary_patterns)
        
    return tuple_to_return


read_data_from_file = open("dataset_30282_4.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

direct_pattern = str(read_strings_from_file[0]).strip()

allowed_distance_input = int(str(read_strings_from_file[1]).strip())

pattern_neighbor_tuple = create_neighbors(direct_pattern, allowed_distance_input,True)

output_file = open("output_possible_neighbor_patterns_list.txt", "w")

output_file.write("direct_pattern_neighbor_list\n\n" + str(pattern_neighbor_tuple[0]).replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", "") + "\n\n")

output_file.write("complementary_pattern_neighbor\n\n" + str(pattern_neighbor_tuple[1]).replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

