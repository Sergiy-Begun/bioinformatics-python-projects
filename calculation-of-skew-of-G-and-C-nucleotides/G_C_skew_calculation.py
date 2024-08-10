#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 13:00:13 2024

@author: sergiybegun
"""

# calculating the difference between occurence of G nucleotide 
# and occurence of C nucleotide 
# to determine the DNA strand and the approximate position of ori

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

investigated_genome_skew = skew_of_G_and_C_occurence("GAGCCACCGCGATA")

print("investigated_genome_skew = ", str(investigated_genome_skew).replace(",","").replace("[", "").replace("]", ""))