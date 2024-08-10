#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 15:37:29 2024

@author: sergiybegun
"""

# calculate the Hamming distance between two genome sequence of equal length

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


read_data_from_file = open("dataset_30278_3.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom_fragments_pair_Hamming_distance = Hamming_distance_calculation(read_strings_from_file[0].strip(), read_strings_from_file[1].strip())

read_data_from_file.close()
    
print("investigated_genom_fragments_pair_Hamming_distance = ", investigated_genom_fragments_pair_Hamming_distance)
    