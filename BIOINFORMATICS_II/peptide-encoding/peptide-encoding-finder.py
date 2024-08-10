#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 08:20:17 2024

@author: sergiybegun
"""

# find all the variants of a given peptide encoding in DNA including reverse compliment


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
        #if (current_amino_acid_element == "STOP"):
        #    break
        amino_acids_string += current_amino_acid_element
    
    return amino_acids_string


def dna_to_rna_translation(dna_string: str):
    """
    translate DNA strand to RNA strand

    Parameters
    ----------
    dna_string : str
        A DNA strand for translation.

    Returns
    -------
    A RNA strand translated from the given DNA strand.

    """
    
    rna_string_output = ""
    #print("len(dna_string) = ", len(dna_string))
    for i in range(0,len(dna_string)):
        current_nucleotide = dna_string[i]
        #print("current_nucleotide = ", current_nucleotide)
        if (current_nucleotide == "t") or (current_nucleotide == "T"):
            rna_string_output += "U"
        else:
            rna_string_output += current_nucleotide
    
    return rna_string_output


def peptide_encodings_finder_function(dna_string: str, peptide_pattern_to_find: str):
    """
    find all variants of peptide encodings in the given dna_string

    Parameters
    ----------
    dna_string : str
        A DNA string for the investigation.
    peptide_pattern_to_fin : str
        A peptide pattern to search in DNA string for the investigation.

    Returns
    -------
    A list of codons including ones from reverse complement strand of DNA.

    """
    
    list_of_codons = []
    
    reverse_compliment_dna_string = reverse_complement(dna_string)[::-1]
    #print("reverse_compliment_dna_string = ", reverse_compliment_dna_string)
    
    length_of_dna_string = len(dna_string)
    
    length_of_codons_pattern_to_find = int(len(peptide_pattern_to_find) * 3)
    
    if (length_of_codons_pattern_to_find > length_of_dna_string):
        raise ValueError("Peptide pattern is too long")
    
    for i in range(0,(length_of_dna_string - length_of_codons_pattern_to_find + 1)):
        print("i = ", i, " of ", length_of_dna_string - length_of_codons_pattern_to_find + 1)
        current_dna_string_fragment = dna_string[i:(i + length_of_codons_pattern_to_find)]
        #print("current_dna_string_fragment = ", current_dna_string_fragment)
        current_reverse_compliment_dna_string_fragment = (reverse_compliment_dna_string[i:(i + length_of_codons_pattern_to_find)])[::-1]
        #print("current_reverse_compliment_dna_string_fragment = ", current_reverse_compliment_dna_string_fragment)
        
        current_rna_string_to_search_codons = dna_to_rna_translation(current_dna_string_fragment)
        current_reverse_compliment_rna_string_to_search_codons = dna_to_rna_translation(current_reverse_compliment_dna_string_fragment)
        
        current_codons_string = rna_translation_function(current_rna_string_to_search_codons)
        
        current_reverse_compliment_codons_string = rna_translation_function(current_reverse_compliment_rna_string_to_search_codons)
        
        
        if (current_codons_string == peptide_pattern_to_find):
            #print("DIRECT")
            #print("current_dna_string_fragment = ", current_dna_string_fragment)
            #print("current_codons_string = ", current_codons_string)
            #print("current_reverse_compliment_dna_string_fragment = ", current_reverse_compliment_dna_string_fragment)
            #print("current_reverse_compliment_codons_string = ", current_reverse_compliment_codons_string)
            
            list_of_codons.append(current_dna_string_fragment)
                    
        if (current_reverse_compliment_codons_string == peptide_pattern_to_find):
            #print("REVERSE")
            #print("current_dna_string_fragment = ", current_dna_string_fragment)
            #print("current_codons_string = ", current_codons_string)
            #print("current_reverse_compliment_dna_string_fragment = ", current_reverse_compliment_dna_string_fragment)
            #print("current_reverse_compliment_codons_string = ", current_reverse_compliment_codons_string)
            
            list_of_codons.append(current_dna_string_fragment)

    
    return list_of_codons





read_data_from_file = open("input_1.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

input_dna_string = str(read_strings_from_file[0]).strip()

peptide_pattern = str(read_strings_from_file[1]).strip()

"""

read_data_from_file = open("Bacillus_brevis.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

input_dna_string = ""

for i in range(0,len(read_strings_from_file)):
    print("reading line ", i, "of ", len(read_strings_from_file))
    input_dna_string += read_strings_from_file[i].strip()

#print("input_dna_string = ", input_dna_string)

peptide_pattern = "VKLFPWFNQY"
print("peptide_pattern = ", peptide_pattern)

"""

list_of_codons_output = peptide_encodings_finder_function(input_dna_string,peptide_pattern)

#print("list_of_codons_output = ", list_of_codons_output)
print("len(list_of_codons_output) = ", len(list_of_codons_output))

output_file = open("output_list_of_codons.txt", "w")

for i in range(0,len(list_of_codons_output)):
    output_file.write((str(list_of_codons_output[i]).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", "") + "\n"))

output_file.close()

read_data_from_file.close()

