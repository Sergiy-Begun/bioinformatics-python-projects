#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:47:25 2024

@author: sergiybegun
"""

def score_motifs(list_of_motifs):
    
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


def most_probable_k_mer_finder(investigated_genom_string: str, k_mer_length: int, probability_matrix):
    """
    

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
    
    if (max_probability == 0.0):
        return []
    
    for current_key in dictionary_of_k_mers.keys():
        if (dictionary_of_k_mers[current_key] == max_probability):
            most_probable_k_mer_list.append(str(current_key))
    
    return most_probable_k_mer_list


def probability_matrix_finder(investigated_genom_list, k_mer_length: int, implement_pseudocounts):
    
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


read_data_from_file = open("input_4.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom_input_list = read_strings_from_file[1].split()

for list_item in investigated_genom_input_list:
    list_item = list_item.strip().capitalize()

k_mer_length_input =  int(str(read_strings_from_file[0]).split()[0].strip())

best_motifs = []
for i in range(0,len(investigated_genom_input_list)):
    best_motifs.append(investigated_genom_input_list[i][0:k_mer_length_input])

best_score = score_motifs(best_motifs)

print("best motifs first version = ", str(best_motifs))

print("best_score = ", best_score)
for k in range(0,(len(investigated_genom_input_list[0]) - k_mer_length_input + 1)):
    cur_first_string_k_mer = investigated_genom_input_list[0][k:(k + k_mer_length_input)]
    print("investigated_genom_input_list[0] = ", investigated_genom_input_list[0])
    print("cur_first_string_k_mer", cur_first_string_k_mer)
    greedy_k_mers = [cur_first_string_k_mer]
    probability_matrix_calculated = probability_matrix_finder(greedy_k_mers,k_mer_length_input,False)


    for i in range (1,len(investigated_genom_input_list)):
        
        print("greedy_k_mers before = ", str(greedy_k_mers))
        print("probability_matrix_calculated before = ", str(probability_matrix_calculated))
        print("i = ",i)
        probability_matrix_calculated = probability_matrix_finder(greedy_k_mers,k_mer_length_input,False)
        print("probability_matrix_calculated after = ", str(probability_matrix_calculated))
        cur_most_probable_k_mer_finder = most_probable_k_mer_finder(investigated_genom_input_list[i], k_mer_length_input, probability_matrix_calculated)
        print("cur_most_probable_k_mer_finder before = ", cur_most_probable_k_mer_finder)
        if (len(cur_most_probable_k_mer_finder) == 0):
            greedy_k_mers.append(investigated_genom_input_list[i][0:k_mer_length_input])
        else:
            greedy_k_mers.append(cur_most_probable_k_mer_finder[0])
            
        print("greedy_k_mers after = ", str(greedy_k_mers))
    cur_score = score_motifs(greedy_k_mers)
    print("cur_score = ", cur_score)
    if (cur_score < best_score):
        best_score = cur_score
        best_motifs = greedy_k_mers
    
    print("best_motifs = ", str(best_motifs))
        

output_file = open("greedy_k_mers_list.txt", "w")

output_file.write(str(best_motifs).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()
