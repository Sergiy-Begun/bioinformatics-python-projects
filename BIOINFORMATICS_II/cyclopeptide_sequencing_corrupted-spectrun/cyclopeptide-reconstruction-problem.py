#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 08:49:20 2024

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
        current_peptide_score = peptide_score_function(current_peptide,spectrum_input,True)
        peptide_score_table.append((current_peptide_score,current_peptide))

    # print("peptide_score_table before = ", sorted(peptide_score_table,reverse=True))

    peptide_score_table = sorted(peptide_score_table,reverse=True)
    
    # print("peptide_score_table after = ", peptide_score_table)
    
    list_of_top_peptides = []
    n_count = 0
    list_filled = False
    
    # print("len(peptide_score_table) = ", len(peptide_score_table))
    """
    for i in range(len(peptide_score_table)):
        current_peptide_score = peptide_score_table[i][0]
        if (list_filled == False):
            list_of_top_peptides.append(peptide_score_table[i][1])
            n_count += 1
            # print("i = ", i)
            # print("n_count = ", n_count, " of ", N_cut)
            # print("current_peptide_score = ", current_peptide_score)
            if (i < (len(peptide_score_table) - 1)):
                next_peptide_score = peptide_score_table[i + 1][0]
                if (n_count >= N_cut) and ((next_peptide_score < current_peptide_score) or (next_peptide_score == 1)):
                    list_filled = True
                    break
            else:
                break
    
    return list_of_top_peptides[0]
"""
    return peptide_score_table[0][1]


def spectrum_formation_for_cyclic_peptide_function(cyclic_peptide_string, cyclic_not_linear: bool, mass_representation: bool):
    """
    give the list of all possible fragment masses including repetitions of the same values

    Parameters
    ----------
    cyclic_peptide_string : str or a list
        A cyclic peptide string for investigation.
    
    cyclic_not_linear : bool
        Is False for cyclic peptide mass spectrometry and True for linear peptide.
        
    mass_representation: bool
        Is True for mass representation of peptide and False for letter representation.

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
        if (mass_representation == False):
            current_element_mass = amino_acid_masses[list_of_starting_variants[0][t]]
        else:
            current_element_mass = list_of_starting_variants[0][t]
        # writing singlets
        list_of_masses.append(current_element_mass)
        full_mass_of_peptide += current_element_mass
    
    #print("list_of_starting_variants = ", list_of_starting_variants)
    
    list_of_masses.append(full_mass_of_peptide)
    
    dictionary_of_masses = {}
    # begin with duplets
    i = 2
    while (i < length_of_cyclic_peptide):
        
        for current_peptide_variant in list_of_starting_variants:
            
            if (mass_representation == False):
            
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
                                
            else:
                current_peptide_string_equivalent = str(current_peptide_variant).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","-").replace("[", "").replace("]", "")

                for j in range(0,(length_of_cyclic_peptide - i + 1)):
                    
                    current_element_for_mass_determination = current_peptide_variant[j:(j + i)]
                    
                    current_element_for_mass_determination_string_equivalent = str(current_element_for_mass_determination).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","-").replace("[", "").replace("]", "")
    
                    count_repetitive_fragments = current_peptide_string_equivalent.count(current_element_for_mass_determination_string_equivalent)
    
                    mass_of_current_element = 0
                    for t in range(0,len(current_element_for_mass_determination)):
                        mass_of_current_element += current_element_for_mass_determination[t]
                                    
                    for repetitive_i in range(0,count_repetitive_fragments):
                        current_repetitive_variant_string_equivalent = current_element_for_mass_determination_string_equivalent + str(repetitive_i)
    
                        if (current_repetitive_variant_string_equivalent not in dictionary_of_masses.keys()):
                            dictionary_of_masses[current_repetitive_variant_string_equivalent] = mass_of_current_element
                        
        i += 1
    
    list_of_masses.extend(list(dictionary_of_masses.values()))
    
    return sorted(list_of_masses)


def peptide_score_function(peptide_string, experimental_spectrum_input: list, mass_representation: bool):
    """
    calculate the score value for a given peptide by comparing
    the experimental and theoretical spectrum and compute the coincidence
    as a score value

    Parameters
    ----------
    peptide_string : str or a list
        A cyclic peptide string for investigation.
    experimental_spectrum_input : list
        The experimental spectrum.
    mass_representation: bool
        Is True for mass representation of peptide and False for letter representation.

    Returns
    -------
    A score value for a given peptide.

    """
    
    peptide_score_value = 0
    
    theoretical_spectrum = sorted(spectrum_formation_for_cyclic_peptide_function(peptide_string,False,True))
    
    # print("theoretical_spectrum = ", theoretical_spectrum)
    length_of_theoretical_spectrum = len(theoretical_spectrum)
    
    length_of_experimental_spectrum = len(experimental_spectrum_input)
    
    # print("length_of_theoretical_spectrum = ", length_of_theoretical_spectrum)
    # print("length_of_experimental_spectrum = ", length_of_experimental_spectrum)
    
    # print("experimental_spectrum_input", experimental_spectrum_input)
    
    already_included = []
    
    for i in range(0,length_of_experimental_spectrum):
        current_experimental_spectrum_element = experimental_spectrum_input[i]
        current_multiplicity = min(theoretical_spectrum.count(current_experimental_spectrum_element),experimental_spectrum_input.count(current_experimental_spectrum_element))
        if (current_experimental_spectrum_element in theoretical_spectrum) and (current_experimental_spectrum_element not in already_included):
            # if current_experimental_spectrum_element < 200:
            #     print("current_multiplicity = ", current_multiplicity, "for the element ", current_experimental_spectrum_element)
            peptide_score_value += (1 * current_multiplicity)
            already_included.append(current_experimental_spectrum_element)
                
    return peptide_score_value


def error_compensation_for_spectrum(current_spectrum_variant: list, max_length_of_list: int, real_spectrum: list):
    """
    reconstruct cyclic peptide by mass spectrum (ONLY for mass representation)

    Parameters
    ----------
    current_spectrum_variant: list
        A zero variant of the list without phantom peaks guaranted.
    
    max_length_of_list : int
        A maximum number of cyclic peptide candidates to remain in the top list.
        
    real_spectrum : list
        peaks in the spectrum with errors.
            
    Returns
    -------
    A top candidate with best score based on list of variants of cyclic peptide.

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
    
    top_candidate = []
    
    current_i = len(current_spectrum_variant[0])
    
    candidate_spectrum = [sorted(copy.deepcopy(list(amino_acid_masses_for_unknown.values())))]
    
    print("candidate_spectrum before = ", str(candidate_spectrum))
    
    for i in range(0,len(candidate_spectrum[0])):
        cur_element_zero_element = candidate_spectrum[0][i]
        # convert the elements of the list into tuples for further using in itertools conversions
        candidate_spectrum[0][i] = (cur_element_zero_element,)
    
    print("candidate_spectrum after = ", str(candidate_spectrum))
    
    for i in range(1,current_i):
        candidate_spectrum.append(current_spectrum_variant)
    
    spectrum_max = max(real_spectrum)
    
    min_amino = min(amino_acid_masses_for_unknown.values()) 
    
    max_i = 1 + int(float(spectrum_max) / float(min_amino))
    
    print("max_i = ", max_i)
    
    for i in range((current_i - 1),max_i):
        last = i
        if len(candidate_spectrum[(i - 1)]) == 0:
            print("I was near last")
            last -= 2
            break
        print("i = ", i)
        #print(time.time() - t1)
        
        copy_of_candidate_spectrum = copy.deepcopy(candidate_spectrum[i])
        
        #print("len(candidate_spectrum[",i,"]) beginning = ", len(candidate_spectrum[i]))
        
        for current_spectrum_element in candidate_spectrum[i]:
            current_total_mass = 0
            #print("current_spectrum_element = ", current_spectrum_element)

            for symb_i in range(0,len(current_spectrum_element)):
                # print("current_spectrum_element[symb_i] = ", current_spectrum_element[symb_i])
                # print("spectrum_max = ", spectrum_max)
                
                current_total_mass += current_spectrum_element[symb_i]
                if (current_total_mass > spectrum_max):
                    copy_of_candidate_spectrum.remove(current_spectrum_element)
                    break
                
            if (current_total_mass == spectrum_max) and (len(top_candidate) == 0):
                top_candidate = current_spectrum_element
                top_candidate_score = peptide_score_function(current_spectrum_element,real_spectrum,True)
                print("top_candidate_score = ", top_candidate_score)
                
            if (current_total_mass == spectrum_max) and (len(top_candidate) > 0):
                current_spectrum_element_score = peptide_score_function(current_spectrum_element,real_spectrum,True)
                if (current_spectrum_element_score > top_candidate_score):
                    top_candidate = current_spectrum_element
                    top_candidate_score = current_spectrum_element_score
                    print("top_candidate_score updated = ", top_candidate_score)
                    
        candidate_spectrum[i] = copy.deepcopy(copy_of_candidate_spectrum)
        #print("len(candidate_spectrum[",i,"]) after cleaning by spectrum comparison = ", len(candidate_spectrum[i]))
        # print("candidate_spectrum[",i,"]) after cleaning by spectrum comparison = ", candidate_spectrum[i])

        candidate_spectrum[i] = copy.deepcopy(sorted(copy_of_candidate_spectrum))
        
        #if (i < 3):
            #print("candidate_spectrum[",i,"]) after duplicate cleaning = ", candidate_spectrum[i])
        # else:
        #     return
        #print("len(candidate_spectrum[",i,"]) after duplicate cleaning = ", len(candidate_spectrum[i]))
        
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
        
        #print("len(candidate_spectrum[",(i + 1),"]) starting values = ", len(candidate_spectrum[i + 1]))
                
        #print(time.time() - t1)
    
    #print("last = ", last)
    
    #print("sorted(list(candidate_spectrum[last])) = ", sorted(list(candidate_spectrum[last])))
    
    # final_variant = copy.deepcopy(trim_peptide_score_table(sorted(list(candidate_spectrum[last])),real_spectrum,max_length_of_list))
    
    return top_candidate


read_data_from_file = open("Tyrocidine_B1_Spectrum_25.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

N_CUT_INPUT = int(str(read_strings_from_file[0]).strip())

spectrum_list = str(read_strings_from_file[1]).split()

for m in range(0,len(spectrum_list)):
    spectrum_list[m] = int(str(spectrum_list[m]).strip())
spectrum_list.remove(0)

#print("spectrum_list = ", spectrum_list)

peptide_variants_output = peptide_reconstruction_from_spectrum(spectrum_list)

print("peptide_variants_output = ", peptide_variants_output)

print("len(peptide_variants_output) = ", len(peptide_variants_output))


cyclic_peptide_top = trim_peptide_score_table(peptide_variants_output, spectrum_list, N_CUT_INPUT)

cyclic_peptide_top = str(cyclic_peptide_top).replace(", ","-").replace("(","").replace(")", "")

print("cyclic_peptide_top after = ", cyclic_peptide_top)

testing_full_spectrum_result = error_compensation_for_spectrum(peptide_variants_output, N_CUT_INPUT, spectrum_list)

print("testing_full_spectrum_result = ", testing_full_spectrum_result)

testing_full_spectrum_result = str(testing_full_spectrum_result).replace(", ","-").replace("(","").replace(")", "")

output_file = open("output_cyclic_peptide_top.txt", "w")

output_file.write(str(testing_full_spectrum_result).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

print(time.time() - t1)
