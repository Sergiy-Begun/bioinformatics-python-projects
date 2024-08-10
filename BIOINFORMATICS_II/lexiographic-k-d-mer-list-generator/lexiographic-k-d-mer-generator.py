#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:50:55 2024

@author: sergiybegun
"""

# generating the list of (k,d)-mers in a lexiographic order from the genom strand

def lexiographic_k_d_mer_generator_function(input_genom_strand: str, k_length: int, d_size: int):
    """
    
    generating the list of (k_length,d_size)-mers in a lexiographic order from the genom strand input_genom_strand

    Parameters
    ----------
    input_genom_strand : str
        A genom strand for cutting.
    k_length : int
        The value of the size of known k-mers at the ends of each cut.
    d_size : int
        The distance between known k-mers at the ends of each cut of the (k,d)-mer.

    Returns
    -------
    A list of (k_length,d_size)-mers of the genom strand input_genom_strand.

    """
    
    list_of_k_d_mers = []
    
    length_of_genom_strand = len(input_genom_strand)
    
    for i in range(0,(length_of_genom_strand - ((2 * k_length) + d_size) + 1)):
        begin_of_the_k_d_mer = input_genom_strand[i:(i + k_length)]
        end_of_the_k_d_mer = input_genom_strand[(i + k_length + d_size):(i + (2 * k_length) + d_size)]
        list_of_k_d_mers.append(str("(" + str(begin_of_the_k_d_mer) + "|" + str(end_of_the_k_d_mer) + ")"))
        #print("list_of_k_d_mers = ", sorted(list_of_k_d_mers))
    
    return sorted(list_of_k_d_mers)

k_d_mers_list_result = lexiographic_k_d_mer_generator_function("GACACATCTCTCA",4,2)

print("k_d_mers_list_result = ", k_d_mers_list_result)

output_file = open("output_k_d_mers_list.txt", "w")

output_file.write(str(k_d_mers_list_result).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()