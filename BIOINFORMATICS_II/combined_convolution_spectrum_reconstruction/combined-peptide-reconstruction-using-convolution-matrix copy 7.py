#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 21:33:12 2024

@author: sergiybegun
"""

import itertools
import copy
import time

t1 = time.time()


def peptide_reconstruction_from_spectrum(allowed_dictionary: list, spectrum: list):
    """
    reconstruct cyclic peptide by mass spectrum

    Parameters
    ----------
    allowed_dictionary : list
        Allowed peptides masses based on convolution.
    spectrum : list
        peaks in the spectrum.

    Returns
    -------
    A list of variants of cyclic peptide.

    """
    
    A_total_mass_was_here = False

    amino_acid_masses_for_unknown = {}
    for i in range(57,201):
        amino_acid_masses_for_unknown[chr(i)] = i
    
    candidate_spectrum = [sorted(copy.deepcopy(allowed_dictionary))]
    
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
        
        #print("len(candidate_spectrum[",i,"]) beginning = ", len(candidate_spectrum[i]))
        
        for current_spectrum_element in candidate_spectrum[i]:
            current_total_mass = 0
            #print("current_spectrum_element = ", current_spectrum_element)
            for symb_i in range(0,len(current_spectrum_element)):
                current_total_mass += current_spectrum_element[symb_i]
                if (current_total_mass not in spectrum):
                    copy_of_candidate_spectrum.remove(current_spectrum_element)
                    break
                
        if (current_total_mass == spectrum_max):
            A_total_mass_was_here = True

        candidate_spectrum[i] = copy.deepcopy(copy_of_candidate_spectrum)
        print("len(candidate_spectrum[",i,"]) after cleaning by spectrum comparison = ", len(candidate_spectrum[i]))

        candidate_spectrum[i] = copy.deepcopy(sorted(copy_of_candidate_spectrum))
        
        #if (i < 3):
            #print("candidate_spectrum[",i,"]) after duplicate cleaning = ", candidate_spectrum[i])
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
    
    if (A_total_mass_was_here == True):
        final_variant = sorted(list(candidate_spectrum[last]))
    else:
        if ((last - 1) >= 1):
            final_variant = sorted(list(candidate_spectrum[last - 1]))
    
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

    amino_acid_masses = {}
    for i in range(57,201):
        amino_acid_masses[chr(i)] = i
    
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

                for j in range(0,(length_of_cyclic_peptide - i + 1)):
                    
                    current_element_for_mass_determination = current_peptide_variant[j:(j + i)]
    
                    count_repetitive_fragments = (str(current_peptide_variant).replace(")", "").replace("(", "").replace(",", "")).count((str(current_element_for_mass_determination)).replace(")", "").replace("(", "").replace(",", ""))
                    #print("count_repetitive_fragments = ", count_repetitive_fragments)
    
                    mass_of_current_element = 0
                    for t in range(0,len(current_element_for_mass_determination)):
                        mass_of_current_element += current_element_for_mass_determination[t]
                                    
                    for repetitive_i in range(0,count_repetitive_fragments):
                        current_repetitive_variant = str(current_element_for_mass_determination) + str(repetitive_i)
                        #print("current_repetitive_variant = ", current_repetitive_variant)
    
                        if (current_repetitive_variant not in dictionary_of_masses.keys()):
                            dictionary_of_masses[current_repetitive_variant] = mass_of_current_element
                        
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
    
    #print("peptide_string = ", peptide_string)
    
    #print("theoretical_spectrum = ", theoretical_spectrum)
    length_of_theoretical_spectrum = len(theoretical_spectrum)
    
    length_of_experimental_spectrum = len(experimental_spectrum_input)
    
    #print("length_of_theoretical_spectrum = ", length_of_theoretical_spectrum)
    #print("length_of_experimental_spectrum = ", length_of_experimental_spectrum)
    
    #print("experimental_spectrum_input", experimental_spectrum_input)
    
    already_included = []
    
    for i in range(0,length_of_experimental_spectrum):
        current_experimental_spectrum_element = experimental_spectrum_input[i]
        current_multiplicity = min(theoretical_spectrum.count(current_experimental_spectrum_element),experimental_spectrum_input.count(current_experimental_spectrum_element))
        #current_multiplicity = experimental_spectrum_input.count(current_experimental_spectrum_element)
        if (current_experimental_spectrum_element in theoretical_spectrum) and (current_experimental_spectrum_element not in already_included):
            # if current_experimental_spectrum_element < 200:
            #     print("current_multiplicity = ", current_multiplicity, "for the element ", current_experimental_spectrum_element)
            peptide_score_value += (1 * current_multiplicity)
            already_included.append(current_experimental_spectrum_element)
                
    return peptide_score_value


def error_compensation_for_spectrum_convolution_based(allowed_dictionary : list, current_spectrum_variant: list, max_length_of_list: int, real_spectrum: list):
    """
    reconstruct cyclic peptide by mass spectrum (ONLY for mass representation)

    Parameters
    ----------
    allowed_dictionary : list
        Allowed peptides masses based on convolution.
    
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

    amino_acid_masses_for_unknown = {}
    for i in range(57,201):
        amino_acid_masses_for_unknown[chr(i)] = i
    
    top_candidate = []
    
    current_i = len(current_spectrum_variant[0])
    
    candidate_spectrum = [sorted(copy.deepcopy(allowed_dictionary))]
    
    #print("candidate_spectrum before = ", str(candidate_spectrum))
    
    for i in range(0,len(candidate_spectrum[0])):
        cur_element_zero_element = candidate_spectrum[0][i]
        # convert the elements of the list into tuples for further using in itertools conversions
        candidate_spectrum[0][i] = (cur_element_zero_element,)
    
    #print("candidate_spectrum after = ", str(candidate_spectrum))
    
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
                
            if (current_total_mass == spectrum_max):
                top_candidate_score = peptide_score_function(current_spectrum_element,real_spectrum,True)
                #print("top_candidate_score = ", top_candidate_score)
                top_candidate.append((top_candidate_score,current_spectrum_element))
            """    
            if (current_total_mass == spectrum_max) and (len(top_candidate) > 0):
                
                
                if (current_spectrum_element_score > top_candidate_score):
                    top_candidate = current_spectrum_element
                    top_candidate_score = current_spectrum_element_score
                    print("top_candidate_score updated = ", top_candidate_score)
            """
                    
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
    
    top_candidate = copy.deepcopy(sorted(top_candidate,reverse=True))
    
    return top_candidate


def cyclopeptide_reconstruction_by_convolution_matrix(input_spectrum: list, most_frequent_max: int):
    """
    
    Parameters
    ----------
    input_spectrum : list
        An experimental spectrum.
        
    most_frequent_max : int
        A number of most frequent peaks to keep in the spectrum as an input for the peptide reconstruction starting point.

    Returns
    -------
    Spectrum of most frequent elements.

    """
    
    output_spectrum = []
    
    spectrum_with_rating = []
    
    full_spectrum = []
    
    convolution_matrix = []
    
    spectrum_length = len(input_spectrum)
    # print("spectrum_length = ", spectrum_length)
    
    for i in range(0,(spectrum_length - 1)):
        convolution_matrix.append([])
        for j in range((i + 1),spectrum_length):
            convolution_matrix[i].append((input_spectrum[j] - input_spectrum[i]))
    
    # print("convolution_matrix = ", convolution_matrix)
    
    for i in range(0,(spectrum_length - 1)):
        for j in range(0,(spectrum_length - i - 1)):
            # print("i = ", i, "j = ", j)
            current_convolution_matrix_element = convolution_matrix[i][j]
            # print("current_convolution_matrix_element[", i," , ", j, "] = ", current_convolution_matrix_element)
            if ((current_convolution_matrix_element >= 57) and (current_convolution_matrix_element <= 200)):
                full_spectrum.append(current_convolution_matrix_element)
    
    full_spectrum_as_a_string = str(full_spectrum)
    for i in range(0, len(full_spectrum)):
        current_sp_element = full_spectrum[i]
        if (current_sp_element not in output_spectrum):
            output_spectrum.append(current_sp_element)
            counts_of_sp_element_in_full_spectrum = full_spectrum_as_a_string.count(str(current_sp_element))
            spectrum_with_rating.append((counts_of_sp_element_in_full_spectrum,current_sp_element))
        
    spectrum_with_rating = sorted(spectrum_with_rating,reverse=True)
    
    #print("spectrum_with_rating = ", spectrum_with_rating)
    
    output_spectrum = []
    spectrum_filled = False
    length_spectrum_with_rating = len(spectrum_with_rating)
    count_of_elements = 0
    while (spectrum_filled == False):
        count_of_elements += 1
        if (count_of_elements <= length_spectrum_with_rating):
            if (count_of_elements <= most_frequent_max):
                output_spectrum.append(spectrum_with_rating[count_of_elements - 1][1])
            else:
                if (spectrum_with_rating[count_of_elements - 2][0] > spectrum_with_rating[count_of_elements - 1][1]):
                    spectrum_filled = True
                else:
                     output_spectrum.append(spectrum_with_rating[count_of_elements - 1][1])   
        else:
            spectrum_filled = True
            break
    
    output_spectrum = copy.deepcopy(sorted(output_spectrum,reverse=False))
    
    #print("output_spectrum = ", output_spectrum)
    
    return output_spectrum


def combined_convolution_spectrum_reconstruction_function(most_frequent_max: int, N_Cut: int, real_spectrum: list):
    """
    Using spectrum comvolution with Leaderboard cyclic peptide reconstruction to get the Top peptide candidate for the spectrum.

    Parameters
    ----------
    most_frequent_max : int
        A number of most frequent peaks to keep in the spectrum as an input for the peptide reconstruction starting point.
    N_Cut : int
        A number of Leaderbord peptides to keep in the peptide reconstruction.
    real_spectrum : list
        A real corrupted spectrum.

    Returns
    -------
    Top candidate peptide in the mass-representation.

    """
    
    top_peptide = ""
    
    reconstructed_variants_from_spectrum = cyclopeptide_reconstruction_by_convolution_matrix(real_spectrum,most_frequent_max)
    
    modified_spectrum_with_refreshed_singlets_part = copy.deepcopy(real_spectrum)
    for i in range(0, len(reconstructed_variants_from_spectrum)):
        current_spect_el = reconstructed_variants_from_spectrum[i]
        if (current_spect_el not in modified_spectrum_with_refreshed_singlets_part):
            modified_spectrum_with_refreshed_singlets_part.append(current_spect_el)
        
    modified_spectrum_with_refreshed_singlets_part = copy.deepcopy(sorted(modified_spectrum_with_refreshed_singlets_part))
    
    # print("modified_spectrum_with_refreshed_singlets_part = ", modified_spectrum_with_refreshed_singlets_part)
    
    top_peptide_candidates = peptide_reconstruction_from_spectrum(reconstructed_variants_from_spectrum,real_spectrum)
    
    # print("top_peptide_candidates = ", top_peptide_candidates)
    
    top_peptide = error_compensation_for_spectrum_convolution_based(reconstructed_variants_from_spectrum, top_peptide_candidates, N_Cut, real_spectrum)
    
    output_top = []
    output_score_max = top_peptide[0][0]
    cur_score = top_peptide[0][0]
    i_count = 0
    while ((cur_score == output_score_max) and (i_count < len(top_peptide))):

        output_top.append(top_peptide[i_count])
        
        i_count += 1
        
        cur_score = top_peptide[i_count][0]
    
    #print("top_peptide = ", top_peptide)
    
    #print("top_peptide[0][1] = ", top_peptide[0][1])
    
    print("output_top = ", output_top)
    
    return output_top

read_data_from_file = open("input_1.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

M_MAX_INPUT = int(str(read_strings_from_file[0]).strip())

N_CUT_INPUT = int(str(read_strings_from_file[1]).strip())

spectrum_list = str(read_strings_from_file[2]).split()

for m in range(0,len(spectrum_list)):
    spectrum_list[m] = int(str(spectrum_list[m]).strip())
#spectrum_list.append(0)

spectrum_list = copy.deepcopy(sorted(spectrum_list))

#print("spectrum_list = ", spectrum_list)

top_candidate = combined_convolution_spectrum_reconstruction_function(M_MAX_INPUT, N_CUT_INPUT, spectrum_list)

output_file = open("output_top_candidate.txt", "w")

for cur_top_candidate in top_candidate:
    
    cur_top_candidate_score = cur_top_candidate[0]
    
    output_file.write(str(cur_top_candidate_score) + "\n")

    cur_top_candidate_peptide = str(cur_top_candidate[1]).replace(", ","-").replace("(","").replace(")", "")
    
    output_file.write(str(cur_top_candidate_peptide).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", "") + "\n")

output_file.close()

read_data_from_file.close()

print(time.time() - t1)
