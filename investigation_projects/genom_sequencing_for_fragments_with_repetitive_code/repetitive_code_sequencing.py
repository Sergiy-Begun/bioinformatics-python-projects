#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 12:20:10 2024

@author: sergiybegun
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# finding positions of DNA fragments with unknown left and/or right ends connections
# for DNA fragments with errors
# based on local minimums of hamming distances

# G_C skew is expected to be constant along strands but should be investigated too


def skew_of_G_and_C_occurence(input_genome_fragment):
    """
    calculation of difference of occurence of G minus C along the genome strand
    from the beginning of the strand to the end

    Parameters
    ----------
    input_genome_fragment : str
        A genome strand sequence fragment for the investigation.

    Returns
    -------
    A list of accumulated differences (G - C) along the strand 
    of length (len(input_genome_fragment) + 1) with first value 0.
    For example, skew_of_G_and_C_occurence("CATGGGCATCGGCCATACGCC") = [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]

    """
    
    list_of_G_C_skew = [0,]
    
    for i in range(0,len(input_genome_fragment)):
        current_genome_nucleotide = input_genome_fragment[i]
        if (current_genome_nucleotide == "C") or (current_genome_nucleotide == "c"):
            list_of_G_C_skew.append((list_of_G_C_skew[i] - 1))
        elif (current_genome_nucleotide == "G") or (current_genome_nucleotide == "g"):
            list_of_G_C_skew.append((list_of_G_C_skew[i] + 1))
        else:
            list_of_G_C_skew.append(list_of_G_C_skew[i])
    
    
    return list_of_G_C_skew


def G_C_skew_minimum_finder(input_genome_strand_fragment):
    """
    finding all the positions of accumulated G-C skew minimum 
    along the genome strand fragment

    Parameters
    ----------
    input_genome_strand_fragment : str
        A genome strand fragment for the investigation.

    Returns
    -------
    A list of all positions along the accumulated G-C skew 
    in the array (list) of length equal to (len(input_genome_strand_fragment) + 1) with first 0 value.
    For example, G_C_skew_minimum_finder("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT") = [11, 24]

    """
    
    list_of_minimums = []
    
    # formation of accumulated skew list with minimum skew determination
    list_of_G_C_skew = [0,]
    
    for i in range(0,len(input_genome_strand_fragment)):
        current_genome_nucleotide = input_genome_strand_fragment[i]
        if (current_genome_nucleotide == "C") or (current_genome_nucleotide == "c"):
            list_of_G_C_skew.append((list_of_G_C_skew[i] - 1))
        elif (current_genome_nucleotide == "G") or (current_genome_nucleotide == "g"):
            list_of_G_C_skew.append((list_of_G_C_skew[i] + 1))
        else:
            list_of_G_C_skew.append(list_of_G_C_skew[i])
        
        if i == 0:
            minimum_skew_value = list_of_G_C_skew[i + 1]
        
        if (i > 0) and (minimum_skew_value > list_of_G_C_skew[i + 1]):
            minimum_skew_value = list_of_G_C_skew[i + 1]
    
    for i in range(1, len(list_of_G_C_skew)):
        if list_of_G_C_skew[i] == minimum_skew_value:
            list_of_minimums.append(i)
    
    return list_of_minimums



investigated_genome_skew = skew_of_G_and_C_occurence("GAGCCACCGCGATA")

print("investigated_genome_skew = ", str(investigated_genome_skew).replace(",","").replace("[", "").replace("]", ""))


read_data_from_file = open("Salmonella_enterica.txt", "r")

read_strings_from_file = read_data_from_file.read().strip().replace("\n","")

G_C_skew_minimum_list = G_C_skew_minimum_finder(read_strings_from_file)

output_file = open("output_minimum_positions_list.txt", "w")

output_file.write(str(G_C_skew_minimum_list).replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

control_output_file = open("control_output_from_reading.txt", "w")

control_output_file.write(read_strings_from_file + "\n" + "length of genome string = " + str(len(read_strings_from_file)))

control_output_file.close()

