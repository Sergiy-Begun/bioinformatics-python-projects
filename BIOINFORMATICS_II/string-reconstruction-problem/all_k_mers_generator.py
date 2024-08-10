#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 18:57:53 2024

@author: sergiybegun
"""

# generate all possible k-mers for a given genom strand

def all_possible_k_mers(genom_strand: str, k_mer_length: int):
    """
    
    give a list of all possible k-mer fragments for a given genom strand

    Parameters
    ----------
    genom_strand : str
        investigated genom strand.
    k_mer_length : int
        k-mer length.

    Returns
    -------
    a list of all possible k-mer fragments for a given genom strand.
    For example, all_possible_k_mers("CAATCCAAC",5) = ["CAATC", "AATCC", "ATCCA", "TCCAA", "CCAAC"]

    """
    
    list_of_k_mers = []
    
    genom_strand_length = len(genom_strand)
    
    if (genom_strand_length < k_mer_length):
        raise ValueError("Genom strand length is shorter than the k-mer length")
    
    for i in range(0,(genom_strand_length - k_mer_length + 1)):
        list_of_k_mers.append(genom_strand[i:(i + k_mer_length)])
    
    return list_of_k_mers

read_data_from_file = open("dataset_30153_3.txt", "r")

read_strings_from_file = read_data_from_file.readlines()

investigated_genom = str(read_strings_from_file[1]).strip()

k_mer_length_input = int(str(read_strings_from_file[0]).strip())

all_possible_k_mers_list = all_possible_k_mers(investigated_genom,k_mer_length_input)

output_file = open("output_all_possible_k_mers_list.txt", "w")

output_file.write(str(all_possible_k_mers_list).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()
