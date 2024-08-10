#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 18:03:33 2024

@author: sergiybegun
"""

# creating the reverse complement DNA sequence

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

read_data_from_file = open("dataset_30273_2.txt", "r")

output_file = open("output_sequence", "w")

read_the_string_from_file = read_data_from_file.read().strip()

reverse_complement_DNA_sequence = reverse_complement(read_the_string_from_file)

output_file.write(reverse_complement_DNA_sequence)

output_file.close()

read_data_from_file.close()