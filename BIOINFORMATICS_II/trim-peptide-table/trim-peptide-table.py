#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 22:56:12 2024

@author: sergiybegun
"""

# compare the experimental and theoretical spectrum and compute the coincidence
# as a score value
# select the N top candidates with the highest score

def trim_peptide_score_table(list_of_peptide_candidates: list, spectrum_input: list, N_cut: int):
    """
    Cut the end of sorted list.

    Parameters
    ----------
    list_of_peptide_candidates : list
        List of peptides to analyze.
    spectrum_input : list
        An experimental spectrum.
    N_cut : int
        A number of top candidates to stay in the list.

    Returns
    -------
    A list of candidates remained.

    """
    
    peptide_score_table = []
    
    length_of_peptide_list = len(list_of_peptide_candidates)
    
    for i in range(0,length_of_peptide_list):
        current_peptide = list_of_peptide_candidates[i]
        current_peptide_score = peptide_score_function(current_peptide,spectrum_input)
        peptide_score_table.append((current_peptide_score,current_peptide))

    print("peptide_score_table before = ", sorted(peptide_score_table,reverse=True))

    peptide_score_table = sorted(peptide_score_table,reverse=True)
    
    print("peptide_score_table after = ", peptide_score_table)
    
    list_of_top_peptides = []
    n_count = 0
    list_filled = False
    for i in range(len(peptide_score_table)):
        current_peptide_score = peptide_score_table[i][0]
        if (list_filled == False):
            list_of_top_peptides.append(str(peptide_score_table[i][1]))
            n_count += 1
            next_peptide_score = peptide_score_table[i + 1][0]
            if (n_count >= N_cut) and (next_peptide_score < current_peptide_score):
                list_filled = True
                break
    
    return list_of_top_peptides


def spectrum_formation_for_cyclic_peptide_function(cyclic_peptide_string: str, cyclic_not_linear: bool):
    """
    give the list of all possible fragment masses including repetitions of the same values

    Parameters
    ----------
    cyclic_peptide_string : str
        A cyclic peptide string for investigation.
    
    cyclic_not_linear : bool
        Is False for cyclic peptide mass spectrometry and True for linear peptide.

    Returns
    -------
    A list of all possible fragment masses including repetitions of the same values
    for a given peptide cyclic_peptide_string.

    """

    amino_acid_masses = {"G": 57,
                        "A": 71,
                        "S": 87,
                        "P": 97,
                        "V": 99,
                        "T": 101,
                        "C": 103,
                        "I": 113,
                        "L": 113,
                        "N": 114,
                        "D": 115,
                        "K": 128,
                        "Q": 128,
                        "E": 129,
                        "M": 131,
                        "H": 137,
                        "F": 147,
                        "R": 156,
                        "Y": 163,
                        "W": 186
                        }
    
    list_of_masses = [0]
    
    list_of_starting_variants = []
    
    # formation of all variants of linear form of cyclic peptide
    
    length_of_cyclic_peptide = len(cyclic_peptide_string)
    if (cyclic_not_linear == False):
        for i in range(0,length_of_cyclic_peptide):
            current_variant = cyclic_peptide_string[i:length_of_cyclic_peptide] + cyclic_peptide_string[0:i]
            list_of_starting_variants.append(current_variant)
    
    if (cyclic_not_linear == True):
        list_of_starting_variants = [cyclic_peptide_string]
    
    full_mass_of_peptide = 0
    for t in range(0,length_of_cyclic_peptide):
        current_element_mass = amino_acid_masses[list_of_starting_variants[0][t]]
        list_of_masses.append(current_element_mass)
        full_mass_of_peptide += current_element_mass
    
    #print("list_of_starting_variants = ", list_of_starting_variants)
    
    list_of_masses.append(full_mass_of_peptide)
    
    dictionary_of_masses = {}
    i = 2
    while (i < length_of_cyclic_peptide):
        
        for current_peptide_variant in list_of_starting_variants:
            for j in range(0,(length_of_cyclic_peptide - i + 1)):
                
                current_element_for_mass_determination = current_peptide_variant[j:(j + i)]

                count_repetitive_fragments = current_peptide_variant.count(current_element_for_mass_determination)

                mass_of_current_element = 0
                for t in range(0,len(current_element_for_mass_determination)):
                    mass_of_current_element += amino_acid_masses[current_element_for_mass_determination[t]]
                                
                for repetitive_i in range(0,count_repetitive_fragments):
                    current_repetitive_variant = current_element_for_mass_determination + str(repetitive_i)

                    if (current_repetitive_variant not in dictionary_of_masses.keys()):
                        dictionary_of_masses[current_repetitive_variant] = mass_of_current_element
                            
                        if mass_of_current_element == 199:
                            print("already there = ", str(current_element_for_mass_determination), "with mass = ", mass_of_current_element)
                    
                    
        i += 1
    
    list_of_masses.extend(list(dictionary_of_masses.values()))
    
    return sorted(list_of_masses)


def peptide_score_function(peptide_string: str, experimental_spectrum_input: list):
    """
    calculate the score value for a given peptide by comparing
    the experimental and theoretical spectrum and compute the coincidence
    as a score value

    Parameters
    ----------
    peptide_string : str
        A cyclic peptide string for investigation.
    experimental_spectrum_input : list
        The experimental spectrum.

    Returns
    -------
    A score value for a given peptide.

    """
    
    peptide_score_value = 0
    
    theoretical_spectrum = sorted(spectrum_formation_for_cyclic_peptide_function(peptide_string,True))
    
    print("theoretical_spectrum = ", theoretical_spectrum)
    length_of_theoretical_spectrum = len(theoretical_spectrum)
    
    length_of_experimental_spectrum = len(experimental_spectrum_input)
    
    print("length_of_theoretical_spectrum = ", length_of_theoretical_spectrum)
    print("length_of_experimental_spectrum = ", length_of_experimental_spectrum)
    
    print("experimental_spectrum_input", experimental_spectrum_input)
    
    already_included = []
    
    for i in range(0,length_of_experimental_spectrum):
        current_experimental_spectrum_element = experimental_spectrum_input[i]
        current_multiplicity = min(theoretical_spectrum.count(current_experimental_spectrum_element),experimental_spectrum_input.count(current_experimental_spectrum_element))
        if (current_experimental_spectrum_element in theoretical_spectrum) and (current_experimental_spectrum_element not in already_included):
            if current_experimental_spectrum_element < 200:
                print("current_multiplicity = ", current_multiplicity, "for the element ", current_experimental_spectrum_element)
            peptide_score_value += (1 * current_multiplicity)
            already_included.append(current_experimental_spectrum_element)
                
    return peptide_score_value

read_data_from_file = open("dataset_30249_3.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

input_peptide_strings = str(read_strings_from_file[0]).split()

list_of_peptide_input = []
for i in range(0,len(input_peptide_strings)):
    list_of_peptide_input.append(str(input_peptide_strings[i]).strip())

experimental_spectrum_from_file = []

for exp_spectrum_el in str(read_strings_from_file[1]).split():
    experimental_spectrum_from_file.append(int(str(exp_spectrum_el).strip()))

experimental_spectrum_from_file = sorted(experimental_spectrum_from_file)

print("input_peptide_strings = ", input_peptide_strings)

print("experimental_spectrum_from_file = ", experimental_spectrum_from_file)

n_input_cutting = int(str(read_strings_from_file[2]).strip())

list_of_remainder = trim_peptide_score_table(list_of_peptide_input, experimental_spectrum_from_file, n_input_cutting)

print("list_of_remainder = ", list_of_remainder)

output_file = open("output_list_of_remainder.txt", "w")

output_file.write(str(list_of_remainder).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()