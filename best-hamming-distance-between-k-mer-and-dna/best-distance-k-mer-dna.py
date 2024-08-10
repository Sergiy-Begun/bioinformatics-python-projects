#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 17:24:07 2024

@author: sergiybegun
"""

# calculation of the best distance between k-mer and DNA strings 
# as a sum of best distances for each string min(k-mer -- current_k-mer_in_current_DNA_string)

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
        raise ValueError("Lengths of genom strands are not equal")
    
    for i in range(0,length_gen_seq_1):
        if (genom_sequence_1[i] != genom_sequence_2[i]):
            Hamming_distance_measure += 1
    
    return Hamming_distance_measure


def best_hamming_distance_k_mer_DNA(k_mer: str, dna_strings):
    """
    finding the minimum hamming distance between investigated k-mer (k_mer)
    and all possible variants of k-mers inside current DNA string
    for each DNA string. The sum of minimum distances is a result. 

    Parameters
    ----------
    k_mer : str
        investigated k-mer.
    dna_strings : list of str
        A list of investigated DNA strings.

    Returns
    -------
    The sum of minimum hamming distances.
    For example, best_hamming_distance_k_mer_DNA("AAA", ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"]) = 5

    """
    
    resulting_hamming_distance = 0
    
    # suppose we have strands of equal lengths 
    # (but it doesn't matter for the rest of the algorithm)
    length_of_genom_strand = len(dna_strings[0])
    
    length_k_mer = len(k_mer)
    
    if (length_of_genom_strand < length_k_mer):
        raise ValueError("Length of genom strand is less than length of k-mer")
    
    for dna_string in dna_strings:
        
        current_min_hamming_distance = (-1)
        
        for i in range(0,(len(dna_string) - length_k_mer + 1)):
            current_dna_k_mer = dna_string[i:(i + length_k_mer)]
            current_hamming_distance = Hamming_distance_calculation(k_mer,current_dna_k_mer)
            
            if (current_min_hamming_distance == (-1)):
                current_min_hamming_distance = current_hamming_distance
            
            if (current_hamming_distance < current_min_hamming_distance):
                current_min_hamming_distance = current_hamming_distance
        
        if (current_min_hamming_distance == (-1)):
            raise ValueError("Check the algorithm")
        
        resulting_hamming_distance += current_min_hamming_distance
    
    return resulting_hamming_distance

read_data_from_file = open("test_input.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom_input_list = read_strings_from_file[1].split()

for list_item in investigated_genom_input_list:
    list_item = list_item.strip().capitalize()

k_mer_input =  read_strings_from_file[0].strip()

best_distance = best_hamming_distance_k_mer_DNA(k_mer_input, investigated_genom_input_list)

print("best distance = ", best_distance)

read_data_from_file.close()
    
    