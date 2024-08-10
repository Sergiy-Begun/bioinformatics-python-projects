#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 13:12:56 2024

@author: sergiybegun
"""
import itertools
import copy
import time

t1 = time.time()

def peptide_reconstruction_from_spectrum(spectrum: list):
    """
    reconstruct cyclic peptide by mass spectrum

    Parameters
    ----------
    spectrum : list
        peaks in the spectrum.

    Returns
    -------
    A list of variants of cyclic peptide.

    """

    amino_acid_masses_for_unknown = {"G": 57,
                                    "A": 71,
                                    "S": 87,
                                    "P": 97,
                                    "V": 99,
                                    "T": 101,
                                    "C": 103,
                                    "I_L": 113,
                                    "N": 114,
                                    "D": 115,
                                    "K_Q": 128,
                                    "E": 129,
                                    "M": 131,
                                    "H": 137,
                                    "F": 147,
                                    "R": 156,
                                    "Y": 163,
                                    "W": 186
                                    }
    
    candidate_spectrum = [sorted(copy.deepcopy(list(amino_acid_masses_for_unknown.values())))]
    
    for i in range(0,len(candidate_spectrum[0])):
        cur_element_zero_element = candidate_spectrum[0][i]
        candidate_spectrum[0][i] = (cur_element_zero_element,)
    
    spectrum_max = max(spectrum)
    
    min_amino = min(amino_acid_masses_for_unknown.values()) 
    
    max_i = 1 + int(float(spectrum_max) / float(min_amino))
    
    print("max_i = ", max_i)

    for i in range(0,max_i):
        last = i
        if len(candidate_spectrum[(i - 1)]) == 0:
            print("I was near last")
            last -= 2
            break
        print("i = ", i)
        print(time.time() - t1)
        
        copy_of_candidate_spectrum = copy.deepcopy(candidate_spectrum[i])
        
        print("len(candidate_spectrum[",i,"]) beginning = ", len(candidate_spectrum[i]))
        
        for current_spectrum_element in candidate_spectrum[i]:
            current_total_mass = 0
            #print("current_spectrum_element = ", current_spectrum_element)
            for symb_i in range(0,len(current_spectrum_element)):
                current_total_mass += current_spectrum_element[symb_i]
                if (current_total_mass not in spectrum):
                    copy_of_candidate_spectrum.remove(current_spectrum_element)
                    break

        candidate_spectrum[i] = copy.deepcopy(copy_of_candidate_spectrum)
        print("len(candidate_spectrum[",i,"]) after cleaning by spectrum comparison = ", len(candidate_spectrum[i]))

        candidate_spectrum[i] = copy.deepcopy(sorted(copy_of_candidate_spectrum))
        
        if (i < 3):
            print("candidate_spectrum[",i,"]) after duplicate cleaning = ", candidate_spectrum[i])
        print("len(candidate_spectrum[",i,"]) after duplicate cleaning = ", len(candidate_spectrum[i]))
        
        candidate_spectrum.append([])
        current_modeled_spectrum_raw = sorted(list(itertools.product(copy.deepcopy(candidate_spectrum[0]),copy.deepcopy(candidate_spectrum[i]),repeat=1)))
        #print("current_modeled_spectrum_raw = ", current_modeled_spectrum_raw)
        for current_modeled_spectrum_el in current_modeled_spectrum_raw:
            current_modeled_spectrum_el = str(current_modeled_spectrum_el).split()
            tupl_conv = []
            for symb_i in range(0,len(current_modeled_spectrum_el)):
                tupl_conv.append(int(current_modeled_spectrum_el[symb_i].replace(")","").replace("(","").replace(",","")))
            candidate_spectrum[i + 1].append(tuple(tupl_conv))
        
        candidate_spectrum[i + 1] = sorted(candidate_spectrum[i + 1])
        
        print("len(candidate_spectrum[",(i + 1),"]) starting values = ", len(candidate_spectrum[i + 1]))
                
        print(time.time() - t1)
    
    print("last = ", last)
    final_variant = sorted(list(candidate_spectrum[last]))    
    
    return final_variant


read_data_from_file = open("dataset_30217_6.txt", "r")

read_strings_from_file = read_data_from_file.read()

spectrum_list = str(read_strings_from_file).split()

for m in range(0,len(spectrum_list)):
    spectrum_list[m] = int(str(spectrum_list[m]).strip())
spectrum_list.remove(0)

#print("spectrum_list = ", spectrum_list)

peptide_variants_output = peptide_reconstruction_from_spectrum(spectrum_list)

#print("peptide_variants_output = ", peptide_variants_output)

print("len(peptide_variants_output) = ", len(peptide_variants_output))

for i in range(0,len(peptide_variants_output)):
    peptide_variants_output[i] = str(peptide_variants_output[i]).replace(", ","-").replace("(","").replace(")", "")

#print("peptide_variants_output after = ", peptide_variants_output)    

output_file = open("output_cyclic_peptide_list.txt", "w")

output_file.write(str(peptide_variants_output).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

print(time.time() - t1)
