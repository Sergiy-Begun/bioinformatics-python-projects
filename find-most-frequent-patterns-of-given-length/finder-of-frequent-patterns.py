#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 14:28:16 2024

@author: sergiybegun
"""

# find most frequent patterns of given length in a string including overlapping of symbols

def max_frequency_finder(investigated_string,seq_length):
    """
    find most frequent patterns of given length in a string including overlapping of symbols

    Parameters
    ----------
    investigated_string : str
        input the investigated string.
    seq_length : int
        the length of the pattern to search in the investigated string.

    Returns
    -------
    the list of most frequent patterns of a given length seq_length in the investigated_string.
    For example max_frequency_finder("ACGTTGCATGTCGCATGATGCATGAGAGCT",4) = ["CATG", "GCAT"]

    """
    full_list_of_patterns = {}
    max_occured_patterns = []
    
    length_of_string = len(investigated_string)
    
    if (seq_length > length_of_string): return []
    
    # formation of the full list of available patterns as a dictionary (keys)
    # and counting their frequency in the investigated_string (dictionary values)
    for i in range(0, (length_of_string - seq_length + 1)):
        current_pattern = investigated_string[i:(i + seq_length)]
        if current_pattern in full_list_of_patterns.keys():
            full_list_of_patterns[current_pattern] += 1
        else:
            full_list_of_patterns[current_pattern] = 1
    
    # finding all the most frequent patterns in the dictionary
    max_count = max(full_list_of_patterns.values())
    
    for current_key in full_list_of_patterns.keys():
        if full_list_of_patterns[current_key] == max_count:
            max_occured_patterns.append(current_key)
    
    return max_occured_patterns


most_frequent_patterns = max_frequency_finder("GATGCCAATATACAAATATACAAACGGTCGCATTCGGTTACGGTCGCGCTCCACCGATATACAAATTCGGTTGATGCCAATATACAAATTCGGTTATTCGGTTACGGTCGCACGGTCGCACGGTCGCGCTCCACCGATTCGGTTACGGTCGCACGGTCGCGCTCCACCGGCTCCACCGGATGCCAGATGCCAACGGTCGCACGGTCGCACGGTCGCGATGCCAACGGTCGCACGGTCGCGCTCCACCGGCTCCACCGACGGTCGCATTCGGTTATATACAAGATGCCAATATACAAGATGCCAGCTCCACCGATATACAAACGGTCGCACGGTCGCATTCGGTTATATACAAACGGTCGCACGGTCGCATATACAAATTCGGTTGATGCCAATTCGGTTGATGCCAGCTCCACCGGCTCCACCGATTCGGTTATATACAAACGGTCGCATTCGGTTGCTCCACCGGCTCCACCGATTCGGTTATTCGGTTATTCGGTTACGGTCGCGATGCCAGCTCCACCGATTCGGTTACGGTCGCGCTCCACCGATTCGGTTACGGTCGCACGGTCGCATTCGGTTACGGTCGCGATGCCAGCTCCACCGGATGCCAGCTCCACCGATTCGGTTATATACAAATATACAAACGGTCGCGATGCCAGCTCCACCGGATGCCAATTCGGTTGATGCCAGCTCCACCGGCTCCACCGACGGTCGCATTCGGTTATATACAAGCTCCACCGATATACAAACGGTCGCACGGTCGCGATGCCAATATACAAGATGCCAACGGTCGCGATGCCAATATACAAATATACAAACGGTCGCGCTCCACCGGATGCCAGCTCCACCGGCTCCACCGACGGTCGCGCTCCACCG",12)

print("most_frequent_patterns = ", most_frequent_patterns)