#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 17:29:01 2024

@author: sergiybegun
"""

# finding the all most frequent k-mers with up to d mismatches in genome strand


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


def most_frequent_k_mer_with_mismatches_finder(input_genome_sequence, k_mer_size, d_allowed_distance, including_reverse_complements):
    """
    finding the all most frequent k-mers of k_mer_size size
    with up to d_allowed_distance mismatches in pattern in input_genome_sequence

    Parameters
    ----------
    input_genome_sequence : str
        Investigated genome sequence.
    k_mer_size : int
        The size of k-mers to find.
    d_allowed_distance : int
        Allowed Hamming distance for such approximate k-mers.
    including_reverse_complements : bool
        With (True) or without (False) finding all possible variants of complementary patterns
        (should both (pattern and complementary) be in the input_genome_sequence).

    Returns
    -------
    The list of all possible most frequent k-mers (k_mer_size) with up to d_allowed_distance mismatches in genome strand
    with (True) or without (False) finding all possible variants of complementary patterns
    (should both (pattern and complementary) be in the input_genome_sequence).
    For example, most_frequent_k_mer_with_mismatches_finder("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1, False) = [ATGC, ATGT, GATG]
    For example, most_frequent_k_mer_with_mismatches_finder("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1, True) = [ATGT, ACAT]

    """
    
    # dictionary for counting real or approximate k-mers individually
    full_dictionary_of_k_mers_counter = {}
        
    length_of_genom = len(input_genome_sequence)
    
    if (length_of_genom < k_mer_size) or (length_of_genom < d_allowed_distance):
        return "Error of input data"
    
    for i in range(0,(length_of_genom - k_mer_size + 1)):
        
        current_sliding_window_fragment = input_genome_sequence[i:(i + k_mer_size)]
        
        current_neighbors_tuple = create_neighbors(current_sliding_window_fragment, d_allowed_distance,False)
        
        current_neighbors_list = current_neighbors_tuple[0]
        current_reverse_complements_neighbors_list = current_neighbors_tuple[1]
            
        for current_direct_pattern in current_neighbors_list:
            
            if current_direct_pattern not in full_dictionary_of_k_mers_counter:
                full_dictionary_of_k_mers_counter[current_direct_pattern] = 1
            else:
                full_dictionary_of_k_mers_counter[current_direct_pattern] += 1
                
        if (including_reverse_complements == True):
            
            for current_complement_pattern in current_reverse_complements_neighbors_list:
                
                if current_complement_pattern not in full_dictionary_of_k_mers_counter:
                    full_dictionary_of_k_mers_counter[current_complement_pattern] = 1
                else:
                    full_dictionary_of_k_mers_counter[current_complement_pattern] += 1

    # Finding maximum counts in k-mers dictionary
    max_count = max(full_dictionary_of_k_mers_counter.values())
    
    list_to_return = []

    for cur_key in full_dictionary_of_k_mers_counter:
        if (full_dictionary_of_k_mers_counter[cur_key] == max_count):
            list_to_return.append(cur_key)
    
    return list_to_return

read_data_from_file = open("dataset_30278_10.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom = str(read_strings_from_file[0]).strip()

string_with_numbers = read_strings_from_file[1].split()

k_mer_input = int(str(string_with_numbers[0]).strip())

allowed_distance = int(str(string_with_numbers[1]).strip())

most_frequent_k_mers_list = most_frequent_k_mer_with_mismatches_finder(investigated_genom, k_mer_input, allowed_distance, True)

output_file = open("output_most_frequent_k_mers_list.txt", "w")

output_file.write(str(most_frequent_k_mers_list).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

# real data computer experiment on Salmonella enterica genome
# we preliminary determined the G-C skew minimum at the positions 3764856 3764858
# will take the average of those close positions (3764857) 
# and will investigate genome strand in +-600 nucleotides to the left, and to the right
read_data_from_file_1 = open("Salmonella_enterica.txt", "r")

read_strings_from_file_1 = read_data_from_file_1.read().strip().replace("\n","")

investigated_genom = str(read_strings_from_file_1).strip()[(3764857 - 600):(3764857 + 600)]

most_frequent_k_mers_list = most_frequent_k_mer_with_mismatches_finder(investigated_genom, 9, 1, True)

output_file_1 = open("output_most_frequent_k_mers_list_Salmonella_enterica.txt", "w")

output_file_1.write(str(most_frequent_k_mers_list).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file_1.close()

read_data_from_file_1.close()

control_output_file = open("control_output_from_reading.txt", "w")

control_output_file.write(read_strings_from_file_1 + "\n" + "length of genome string = " + str(len(read_strings_from_file_1)))

control_output_file.close()

