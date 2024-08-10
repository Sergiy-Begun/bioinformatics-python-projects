#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 21:42:14 2024

@author: sergiybegun
"""

# count all variants of DNA combinations for cyclic peptide


def cyclic_peptide_counts(amino_acids_string_string: str):
    
    rna_translation_dictionary = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
        "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
        "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
        "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
        "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
        "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
        "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
        "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
        "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    count_dictionary =  {"F": 2,
                        "L": 6,
                        "S": 6,
                        "Y": 2,
                        "STOP": 3,
                        "C": 2,
                        "W": 1,
                        "P": 4,
                        "H": 2,
                        "Q": 2,
                        "R": 6,
                        "I": 3,
                        "M": 1,
                        "T": 4,
                        "N": 2,
                        "K": 2,
                        "V": 4,
                        "A": 4,
                        "D": 2,
                        "E": 2,
                        "G": 4}
    """
    count_dictionary = {}
    
    list_of_rna_translation_dictionary_values = list(rna_translation_dictionary.values())
    
    for i in range(0, len(list_of_rna_translation_dictionary_values)):
        current_value = list_of_rna_translation_dictionary_values[i]
        if (current_value in count_dictionary):
            count_dictionary[current_value] += 1
        else:
            count_dictionary[current_value] = 1
    
    print("count_dictionary = ", count_dictionary)
    """
    
    number_of_dna_variants = 1
    
    for i in range(0,len(amino_acids_string_string)):
        number_of_dna_variants *= count_dictionary[amino_acids_string_string[i]]
    
    return number_of_dna_variants



read_data_from_file = open("input_1.txt", "r")

read_strings_from_file = read_data_from_file.read().strip()

input_amino_acids_string = read_strings_from_file

number_of_variants = cyclic_peptide_counts(input_amino_acids_string)

print("number_of_variants = ", number_of_variants)

read_data_from_file.close()
