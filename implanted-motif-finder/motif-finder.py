#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 00:56:49 2024

@author: sergiybegun
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 17:29:01 2024

@author: sergiybegun
"""

# finding the all most frequent k-mers with up to d mismatches in genome strand

import copy

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
        raise ValueError("Error = length mismatch")
    
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
            raise ValueError("Error = " + current_symbol + " i = ", i)
        
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


def motif_finder(list_of_input_genome_sequences, k_mer_size: int, d_allowed_distance: int, including_reverse_complements):
    """
    finding the all mutated variants of motif
    that should appear in ALL the strings with a given level of mutations

    Parameters
    ----------
    list_of_input_genome_sequences : list
        Investigated genome sequences.
    k_mer_size : int
        The size of k-mers to find.
    d_allowed_distance : int
        Allowed Hamming distance for such approximate k-mers.
    including_reverse_complements : bool
        With (True) or without (False) finding all possible variants of complementary patterns
        (should both (pattern and complementary) be in the input_genome_sequence).

    Returns
    -------
    The list of all possible most frequent motif k-mers (k_mer_size) with up to d_allowed_distance mismatches in genome strands
    with (True) or without (False) finding all possible variants of complementary patterns.
    For example, motif_finder("AAAAA AAAAA AAAAA", 3, 1, True) = [AAA, AAC, AAG, AAT, ACA, AGA, ATA, CAA, GAA, TAA]

    """    
    list_of_final_motif_candidates = []
    
    full_set_of_first_string_k_mers = []
    
    first_string_input_genome_sequence = list_of_input_genome_sequences[0].strip()
    length_of_first_string_of_the_genom = len(first_string_input_genome_sequence)
        
    if (length_of_first_string_of_the_genom < k_mer_size) or (length_of_first_string_of_the_genom < d_allowed_distance):
        raise ValueError("Error of input data")
    
    for i in range(0,(length_of_first_string_of_the_genom - k_mer_size + 1)):
        
        current_sliding_window_fragment = first_string_input_genome_sequence[i:(i + k_mer_size)]
        
        current_neighbors_tuple = create_neighbors(current_sliding_window_fragment, d_allowed_distance,False)
        
        current_neighbors_list = current_neighbors_tuple[0]
        
        for first_string_current_pattern in current_neighbors_list:
            if first_string_current_pattern not in full_set_of_first_string_k_mers:
                full_set_of_first_string_k_mers.append(str(first_string_current_pattern))
    
    if (len(list_of_input_genome_sequences) <= 1):
        return []
    
    list_of_final_motif_candidates = copy.deepcopy(full_set_of_first_string_k_mers)
        
    n = 0
    at_least_one_connection = False
    while (n < len(list_of_final_motif_candidates)):
        
        for i in range(1,len(list_of_input_genome_sequences)):
            at_least_one_connection = False
            genom_string = list_of_input_genome_sequences[i]
            length_of_current_genom_string = len(genom_string)
                        
            if (len(list_of_final_motif_candidates) == 0):
                break
    
            for j in range(0,(length_of_current_genom_string - k_mer_size + 1)):

                current_sliding_window_fragment_per_string = genom_string[j:(j + k_mer_size)]
                
                current_distance = Hamming_distance_calculation(current_sliding_window_fragment_per_string,list_of_final_motif_candidates[n])
                
                if (current_distance <= d_allowed_distance):
                    at_least_one_connection = True
        
            if (at_least_one_connection == False):
                list_of_final_motif_candidates.pop(n)
                if (n > (-1)):
                    n -= 1
            
        n += 1  
    
    return list_of_final_motif_candidates

read_data_from_file = open("dataset_30302_8.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom_list = str(read_strings_from_file[1]).split()

string_with_numbers = read_strings_from_file[0].split()

k_mer_input = int(str(string_with_numbers[0]).strip())

allowed_distance = int(str(string_with_numbers[1]).strip())

hidden_message_list = motif_finder(investigated_genom_list, k_mer_input, allowed_distance, True)

output_file = open("hidden_message_list.txt", "w")

output_file.write(str(hidden_message_list).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

# finding minimum

read_data_from_file_1 = open("dataset_30304_9.txt", "r")

read_strings_from_file_1 = read_data_from_file_1.readlines()

k_mer_input_1 = int(str(read_strings_from_file_1[0]).strip())

investigated_genom_list_1 = []
for m in range(1,len(read_strings_from_file_1)):
    investigated_genom_list_1.append(str(read_strings_from_file_1[m]).strip())

allowed_distance_1 = 0

list_of_candidates_for_minimizing = []

while (len(list_of_candidates_for_minimizing) == 0):

    list_of_candidates_for_minimizing = motif_finder(investigated_genom_list_1, k_mer_input_1, allowed_distance_1, True)
    
    allowed_distance_1 += 1

output_file_1 = open("minimum_finder_list.txt", "w")

output_file_1.write(str(list_of_candidates_for_minimizing).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file_1.close()

read_data_from_file_1.close()


