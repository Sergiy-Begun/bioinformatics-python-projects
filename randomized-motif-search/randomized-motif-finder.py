#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 22:57:07 2024

@author: sergiybegun
"""

# find best motifs of the DNA strings by randomly selecting the k-mers in those string
# and scoring of the obtained sets of motif to find the set with the best score

import numpy as np
import copy


def most_probable_k_mer_finder(investigated_genom_string: str, k_mer_length: int, probability_matrix):
    """
    finding most probable k-mer in a given genom string based on probbability matrix (4 x k)

    Parameters
    ----------
    investigated_genom_string : str
        A genom string to find the most probable k-mer inside.
    k_mer_length : int
        A k-mer size.
    probability_matrix : TYPE
        An array of probabilities [4 x k].

    Returns
    -------
    A list of most probable k-mers in a given genom string.

    """
    
    most_probable_k_mer_list = []
    
    dictionary_to_convert_to_nucleotides_probability_matrix_strings = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3
        }
    
    dictionary_of_k_mers = {}
    
    length_of_genom_string = len(investigated_genom_string)
    
    if (k_mer_length > length_of_genom_string):
        raise ValueError("The k-mer length is greater than the genom length")
    
    for i in range(0,(length_of_genom_string - k_mer_length + 1)):
        
        current_moving_window_k_mer = investigated_genom_string[i:(i + k_mer_length)]
        
        current_probability = 1.0
        for j in range(0,k_mer_length):
            decoded_row_number_for_probability_matrix = dictionary_to_convert_to_nucleotides_probability_matrix_strings[current_moving_window_k_mer[j]]
            current_decoded_probability_for_nucleotide = probability_matrix[decoded_row_number_for_probability_matrix][j]
            
            current_probability *= current_decoded_probability_for_nucleotide
        
        if (current_moving_window_k_mer not in dictionary_of_k_mers.keys()):
            dictionary_of_k_mers[current_moving_window_k_mer] = current_probability
    
    max_probability = max(dictionary_of_k_mers.values())
    
    for current_key in dictionary_of_k_mers.keys():
        if (dictionary_of_k_mers[current_key] == max_probability):
            most_probable_k_mer_list.append(str(current_key))
    
    return most_probable_k_mer_list


def probability_matrix_finder(investigated_genom_list, k_mer_length: int, implement_pseudocounts):
    """
    Generate the probability matrix for nucleotides versus nucleotides' position

    Parameters
    ----------
    investigated_genom_list : list of str
        A list of investigated genom strands.
    k_mer_length : int
        The length of k-mer for which the matrix is calculated 
        from zero position of current genon string to (k_mer_length - 1) position.
    implement_pseudocounts : bool
        If True we use non-zero filling value for initialization filling of probability matrix.

    Raises
    ------
    ValueError
        If length of k-mer greater than the length of genom strand.

    Returns
    -------
    probability_matrix : [4 x k] of float
        The probability to find the certain nucleotide (A, C, G, or T) 
        at the certain position of genom string.

    """
    
    dictionary_to_convert_to_nucleotides_probability_matrix_strings = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3
        }
        
    # matrix [4 x k_mer_length]
    probability_matrix = [
        [],
        [],
        [],
        []
        ]
    
    if (implement_pseudocounts == False):
        filling_value = 0.0
        
    if (implement_pseudocounts == True):
        filling_value = 1.0
    
    # initialization of probability matrix
    for i in range(0,4):
        for j in range(0,k_mer_length):
            probability_matrix[i].append(filling_value)
    
    number_of_dna_strings = len(investigated_genom_list)
    
    accumulated_number_of_nucleotides_investigated = []
    for j in range(0,k_mer_length):
        accumulated_number_of_nucleotides_investigated.append(0.0)
    
    # all are of equal length
    length_of_genom_string = len(investigated_genom_list[0])
    
    if (length_of_genom_string < k_mer_length):
        raise ValueError("the length of genome string is less than k-mer length")
    
    for i in range(0,number_of_dna_strings):
        for j in range(0,k_mer_length):
            accumulated_number_of_nucleotides_investigated[j] += 1.0
            converted_to_number_position_current_nucleotide = dictionary_to_convert_to_nucleotides_probability_matrix_strings[investigated_genom_list[i][j]]
            probability_matrix[converted_to_number_position_current_nucleotide][j] += 1.0
    
    for i in range(0,4):
        for j in range(0,k_mer_length):
            probability_matrix[i][j] /= accumulated_number_of_nucleotides_investigated[j]
    
    return probability_matrix


def score_motifs(list_of_motifs):
    """
    Calculate the score of motifs based on the frequency of nucleotide 
    in the columns of matrix constructed by placement as rows the motifs
    from the list_of_motifs.

    Parameters
    ----------
    list_of_motifs : list of str
        The list of k-mers of equal length
        which are the current version of best motifs.

    Returns
    -------
    score_result : int
        The calculated sum of values per column of constructed matrix.

    """
    
    length_of_motif = len(list_of_motifs[0])
    
    number_of_motifs = len(list_of_motifs)
    
    most_frequent_nucleotides = {}
    
    score_result = 0
    
    for i in range(0,length_of_motif):
        dictionary_to_count = {
            "A": 0,
            "C": 0,
            "G": 0,
            "T": 0
            }
        for j in range(0,number_of_motifs):
            dictionary_to_count[list_of_motifs[j][i].capitalize()] += 1
        
        max_count = max(dictionary_to_count.values())
        
        for dic_key in dictionary_to_count.keys():
            if (dictionary_to_count[dic_key] == max_count):
                current_most_frequent_nucleotide = str(dic_key)
        
        most_frequent_nucleotides[i] = (current_most_frequent_nucleotide,max_count)
    
    for i in range(0,length_of_motif):
        score_result += (number_of_motifs - most_frequent_nucleotides[i][1])
    
    return score_result


def randomized_motif_finder(investigated_DNA_strings, k_mer_length: int, number_of_runs: int, Gibbs_sampling):
    """
    
    Finding a set of motifs by using random methods

    Parameters
    ----------
    investigated_DNA_strings : list of str
        The list of investigated DNA strings.
    k_mer_length : int
        The length of k-mer for motifs to search inside investigated DNA strings.
    number_of_runs : int
        The number of randomized algorithm runs.
    Gibbs_sampling: bool
        If True we use Gibbs sampling algorithm, if False we use randomized algorithm.

    Returns
    -------
    The list of best candidates for motifs.

    """
    
    counter = 0
    best_motif_score = 10000000000
    best_motif_list_to_find = []
    
    if Gibbs_sampling == True:
        
        randomly_generated_current_motif_list = []
        
        for i in range(0,len(investigated_DNA_strings)):
            
            current_DNA_string = investigated_DNA_strings[i]
            length_of_current_DNA_string = len(current_DNA_string)
            
            if (length_of_current_DNA_string < k_mer_length):
                raise ValueError("The length of k-mer is greater than the length of genom strand")
                
            rng = np.random.default_rng()
            
            current_random_starting_position = int(float(length_of_current_DNA_string - k_mer_length + 1) * rng.random())
            
            # print("current_random_starting_position = ", current_random_starting_position)
            
            randomly_generated_current_motif_list.append(current_DNA_string[current_random_starting_position:(current_random_starting_position + k_mer_length)])
            
            if (len(randomly_generated_current_motif_list[i]) < k_mer_length):
                raise ValueError("Generated defective k-mer")
        
        current_probability_matrix = probability_matrix_finder(randomly_generated_current_motif_list, k_mer_length, True)
        
        current_motif_list = randomly_generated_current_motif_list
        
        for t in range(0,number_of_runs):
            
            current_motif_list_score = score_motifs(current_motif_list)
            print("current_motif_list before = ", str(current_motif_list))
            random_string_number = int(rng.random() * (len(investigated_DNA_strings) - 1))
            temp_matrix = copy.deepcopy(current_motif_list).pop(random_string_number)
            print("temp_matrix = ", temp_matrix)
            current_DNA_string_most_probable_k_mer = most_probable_k_mer_finder(investigated_DNA_strings[random_string_number], k_mer_length, temp_matrix)[0]
            
            current_motif_list_score = score_motifs(current_motif_list)
            current_motif_list[random_string_number] = current_DNA_string_most_probable_k_mer
            
            print("current_motif_list after = ", str(current_motif_list))
            if (float(best_motif_score) > float(current_motif_list_score)):
                best_motif_list_to_find = current_motif_list
                best_motif_score = current_motif_list_score
            
            print("best_motif_score after = ", best_motif_score)
            
            print("best_motif_list_to_find after = ", best_motif_list_to_find)

    if Gibbs_sampling == False:
        for m in range(0,number_of_runs):
    
            
            current_run_in_percents = (float(m) / float(number_of_runs)) * 100.0
            
            print("current run position = ", current_run_in_percents, " %")
            
            randomly_generated_current_motif_list = []
            
            for i in range(0,len(investigated_DNA_strings)):
                
                current_DNA_string = investigated_DNA_strings[i]
                length_of_current_DNA_string = len(current_DNA_string)
                
                if (length_of_current_DNA_string < k_mer_length):
                    raise ValueError("The length of k-mer is greater than the length of genom strand")
                    
                rng = np.random.default_rng()
                
                current_random_starting_position = int(float(length_of_current_DNA_string - k_mer_length + 1) * rng.random())
                
                # print("current_random_starting_position = ", current_random_starting_position)
                
                randomly_generated_current_motif_list.append(current_DNA_string[current_random_starting_position:(current_random_starting_position + k_mer_length)])
                
                if (len(randomly_generated_current_motif_list[i]) < k_mer_length):
                    raise ValueError("Generated defective k-mer")
            
            current_probability_matrix = probability_matrix_finder(randomly_generated_current_motif_list, k_mer_length, True)
        
            consensus = False
            while consensus == False:
                
                current_motif_list = []
                
                for i in range(0,len(investigated_DNA_strings)):
                    current_DNA_string_most_probable_k_mer = most_probable_k_mer_finder(investigated_DNA_strings[i], k_mer_length, current_probability_matrix)[0]
                    
                    current_motif_list.append(current_DNA_string_most_probable_k_mer)                
            
                print("current_motif_list before = ", str(current_motif_list))
                current_motif_list_score = score_motifs(current_motif_list)
                
                print("current_motif_list_score before = ", current_motif_list_score)
                
                print("best_motif_score before = ", best_motif_score)
                
                print("best_motif_list_to_find before = ", best_motif_list_to_find)
                
                counter += 1
                        
                if (float(best_motif_score) > float(current_motif_list_score)):
                    best_motif_list_to_find = current_motif_list
                    best_motif_score = current_motif_list_score
                else:
                    consensus = True
                
                current_probability_matrix = probability_matrix_finder(current_motif_list, k_mer_length, True)
        
                print("best_motif_score after = ", best_motif_score)
                
                print("best_motif_list_to_find after = ", best_motif_list_to_find)
    
    return best_motif_list_to_find

read_data_from_file = open("dataset_30309_11.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom_input_list = read_strings_from_file[1].split()

for list_item in investigated_genom_input_list:
    list_item = list_item.strip().capitalize()

k_mer_length_input =  int(str(read_strings_from_file[0]).split()[0].strip())

# the number of randomized algorithm runs to find solution
number_of_sampling = 3000

best_motifs_list = randomized_motif_finder(investigated_genom_input_list, k_mer_length_input, number_of_sampling, False)

output_file = open("output_best_motifs_list.txt", "w")

output_file.write(str(best_motifs_list).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()
