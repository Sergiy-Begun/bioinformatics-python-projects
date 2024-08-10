#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 12:36:29 2024

@author: sergiybegun
"""

# count of pattern occurence in the string including overlapping of symbols

def count_overlapped_occurence(investigated_string,pattern_seq):
    """
    count of pattern occurence in the string including overlapping of symbols

    Parameters
    ----------
    investigated_string : str
        input the investigated string.
    pattern_seq : str
        pattern to search in the investigated string.

    Returns
    -------
    the occurence number of pattern_seq in the investigated_string.
    For example count_overlapped_occurence("rimerimer","rimer") = 2

    """
    occurence_counter = 0
    length_of_pattern = len(pattern_seq)
    length_of_string = len(investigated_string)
    
    if (length_of_pattern > length_of_string): return 0
    
    for i in range(0, (length_of_string - length_of_pattern + 1)):
        pattern_finder = investigated_string.find(pattern_seq,i,(i + length_of_pattern))
        if pattern_finder >= 0:
            occurence_counter += 1
    
    return occurence_counter
 
inv_count = count_overlapped_occurence("GGATTACGTTTCTCCTTGAGAGTCGGTCACTTCTCCTTTCTCCTTTCTCCTGTTCTCCTTTCTCCTCTTTCTCCTTTCTCCTATTCTCCTTTCTCCTCATTTGCTTCTCCTCATATTCTCCTGTTCTCCTTTCTCCTTTAGGTTCTGCCGGAGTATTCTCCTCAAGAATTTCTCCTTCTTCTCCTAAGATTTCTCCTTTCTCCTGACTAGCTTCTCCTTTCTCCTCGATTCTCCTGCGCTTCTCCTGGCGTTTCTCCTCATTCTCCTGTTTCTCCTGTCTTTCTCCTTTCTCCTGACGTGGAACTTCTCCTGTTCTCCTTCGATTCTCCTAAGTTCTCCTGCTTCTCCTTTTCTCCTGTGTTTAACGTTCTCCTTTCTCCTCTTACTTCTCCTGTTGCTTCTCCTGTTCTCCTATATTTTCTCCTGCGTTTCTCCTTTCTCCTCTTCTCCTTTTCTCCTATCGACAAATTCTCCTAGGTTCTCCTGATTCTCCTTTCTCCTTCGAGTTTCTCCTTAAGTTCTCCTTTCTCCTCCCTGCATTTCTCCTGTTTACCTTCTCCTCTTTTCTCCTCTACGTTTCTCCTTTCTCCTTACGACTACCTTCTCCTCAATTCTCCTATTCTCCTATTCTCCTGAGATTCTCCTTTCTCCTACGTTCTCCTTTCTCCTTTCTCCTTTCTCCTAGGGTTCTCCTTCGTGTGGTCTTCTCCTTGCATTCTCCTTTCTCCTATTCTCCTTTCTCCTGACCTTCTCCTTTCTCCTTTCTCCTTTCTCCTTTCTCCTGCGCGTTCTCCTGTTCTCCTATGTTCTCCTTTCTCCTCAAATCTTCTCCTACTCGACTTCTCCTTTCTCCTTTCTCCTTTCTCCTTTCTCCTGTCGTCTTCTCCTTTCTCCTTTCTCCTATTCTCCTTTCTCCTATTCTCCTCGAGTAGTTCTCCTTTCGCTTCTCCTCTTTCTCCTCTTTCTCCTCTTCTCCTTCTTGCATTTCTCCTAGCTTCTCCT","TTCTCCTTT")

print("inv_count = ", inv_count)