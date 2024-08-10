#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 21:42:14 2024

@author: sergiybegun
"""

# translate 3-mers in the RNA string into amino acids string using a dictionary of codons
# stop codons should be ignored


def rna_translation_function(rna_string: str):
    """
    translate 3-mers in the RNA string into amino acids

    Parameters
    ----------
    rna_string : str
        RNA string for translation.

    Returns
    -------
    A string of aminoacids coded in the RNA string (rna_string).
    For example, rna_translation_function("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA") = "MAMAPRTEINSTRING"

    """
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
    
    amino_acids_string = ""
    
    for i in range(0,len(rna_string),3):
        current_3_mer = rna_string[i:(i + 3)]
        current_amino_acid_element = rna_translation_dictionary[current_3_mer]
        if (current_amino_acid_element != "STOP"):
            amino_acids_string += current_amino_acid_element    
    
    return amino_acids_string



read_data_from_file = open("dataset_30213_4.txt", "r")

read_strings_from_file = read_data_from_file.read().strip()

input_rna_string = read_strings_from_file

amino_acids_string = rna_translation_function(input_rna_string)

print("amino_acids_string = ", amino_acids_string)

output_file = open("output_amino_acids_string.txt", "w")

output_file.write((str(amino_acids_string).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", "") + "\n"))

output_file.close()

read_data_from_file.close()
