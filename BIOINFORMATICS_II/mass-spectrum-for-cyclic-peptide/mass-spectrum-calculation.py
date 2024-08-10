#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 19:22:30 2024

@author: sergiybegun
"""

# generate all possible combinations of cyclic peptide masses

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
                    
                mass_of_current_element = 0
                for t in range(0,len(current_element_for_mass_determination)):
                    mass_of_current_element += amino_acid_masses[current_element_for_mass_determination[t]]
                                
                if (current_element_for_mass_determination not in dictionary_of_masses.keys()):
                    dictionary_of_masses[current_element_for_mass_determination] = mass_of_current_element
        i += 1
    
    list_of_masses.extend(list(dictionary_of_masses.values()))
    
    return sorted(list_of_masses)



read_data_from_file = open("dataset_30248_2.txt", "r")

read_strings_from_file = read_data_from_file.read()

input_peptide_string = str(read_strings_from_file).strip()

print("input_peptide_string = ", input_peptide_string)

list_of_masses_output = spectrum_formation_for_cyclic_peptide_function(input_peptide_string,True)

print("list_of_masses_output = ", list_of_masses_output)
print("len(list_of_masses_output) = ", len(list_of_masses_output))

output_file = open("output_list_of_masses.txt", "w")

output_file.write((str(list_of_masses_output).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", "") + "\n"))

output_file.close()

read_data_from_file.close()


