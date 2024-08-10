#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 22:13:16 2024

@author: sergiybegun
"""

import time
#import pickle
import concurrent.futures

t1 = time.time()

# count all possible variants of peptides for a given total mass of peptide


def counter_of_peptide_variants(maximum_mass: int):
    """
    count all possible variants of peptides for a given total mass of peptide

    Parameters
    ----------
    maximum_mass : int
        A maximum mass for a resulting peptide.

    Returns
    -------
    A number of possible combinations.

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
    
    relative_variations_masses = [0, 14, 30, 40, 42, 44, 46, 56, 57, 58, 71, 72, 74, 80, 90, 99, 106, 129]
    
    #list_of_amino_acid_masses_for_unknown_keys = list(amino_acid_masses_for_unknown.keys())
    
    given_mass_dictionary = {maximum_mass: 0}
    
    #length_addons = len(amino_acid_masses_for_unknown.keys())
    
    length_relative = len(relative_variations_masses)
    
    thresholds = {}
    
    minimum_addon = min(amino_acid_masses_for_unknown.values())
    
    maximum_addon = max(amino_acid_masses_for_unknown.values())
    
    current_maximum_length_of_peptide = int(float(maximum_mass) / float(minimum_addon)) + 1
    current_minimum_length_of_peptide = int(float(maximum_mass) / float(maximum_addon)) - 1
    
    thresholds = {"minimum_start": current_minimum_length_of_peptide, "maximum_end": current_maximum_length_of_peptide}
    
    print("thresholds = ", thresholds)
            
    print("given_mass_dictionary[maximum_mass] = ", given_mass_dictionary[maximum_mass])
    
    i =  thresholds["minimum_start"]
    
    while (i <= thresholds["maximum_end"]):
        monit = {"counter": 0}
        # formation of peptide dictionary
        print("i = ", i)
        
        starting_variant_nucletide = ""
        for j in range(0,i):
            starting_variant_nucletide += "G"
        
        starting_variant_nucleotide_as_a_dict = {starting_variant_nucletide: [(i * 57)]}
        
        print("starting_variant_nucleotide_as_a_dict = ", starting_variant_nucleotide_as_a_dict)
        
        m0_starting_mass = starting_variant_nucleotide_as_a_dict[starting_variant_nucletide][0]
        
        if (m0_starting_mass == maximum_mass):
            given_mass_dictionary[maximum_mass] += 1
        
        if (m0_starting_mass >= maximum_mass):
            return given_mass_dictionary[maximum_mass]
        
        for j in range(1,length_relative):
            starting_variant_nucleotide_as_a_dict[starting_variant_nucletide].append((m0_starting_mass + relative_variations_masses[j]))
        
        maximum_total = (20 ** i)
        
        def adding_other_positions(m0: int):
            
            if ((monit["counter"] / 100000) == int((monit["counter"] / 100000))):
                print(monit["counter"], "\t", maximum_total)
                
            added_here = 0
            
            for i_1 in range(0,length_relative):
                
                mass_of_new_element = m0
                
                if ((monit["counter"] / 100000) == int((monit["counter"] / 100000))):
                    print(monit["counter"], "\t", maximum_total)
                
                if i <= 1:
                    break
                if (i > 1):
                    delta_1 = relative_variations_masses[i_1]
                    monit["counter"] += 1                    
                    if ((mass_of_new_element + delta_1) == maximum_mass):
                        given_mass_dictionary[maximum_mass] += 1
                        added_here += 1
                    if (mass_of_new_element > maximum_mass):
                        
                        print("break by mass parameter = ", mass_of_new_element)
                        break
                for i_2 in range(0,length_relative):
                    if i <= 2:
                        break
                    if (i > 2):
                        delta_2 = relative_variations_masses[i_2]
                        monit["counter"] += 1
                        if ((mass_of_new_element + delta_1 + delta_2) == maximum_mass):
                            given_mass_dictionary[maximum_mass] += 1
                            added_here += 1
                        if (mass_of_new_element > maximum_mass):
                            
                            print("break by mass parameter = ", mass_of_new_element)
                            break
                    for i_3 in range(0,length_relative):
                        if i <= 3:
                            break
                        if (i > 3):
                            delta_3 = relative_variations_masses[i_3]
                            monit["counter"] += 1
                            if ((mass_of_new_element + delta_1 + delta_2 + delta_3) == maximum_mass):
                                given_mass_dictionary[maximum_mass] += 1
                                added_here += 1
                            if (mass_of_new_element > maximum_mass):
                                
                                print("break by mass parameter = ", mass_of_new_element)
                                break
                        for i_4 in range(0,length_relative):
                            if i <= 4:
                                break
                            if (i > 4):
                                delta_4 = relative_variations_masses[i_4]
                                monit["counter"] += 1
                                if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4) == maximum_mass):
                                    given_mass_dictionary[maximum_mass] += 1
                                    added_here += 1
                                if (mass_of_new_element > maximum_mass):
                                    
                                    print("break by mass parameter = ", mass_of_new_element)
                                    break
                            for i_5 in range(0,length_relative):
                                if i <= 5:
                                    break
                                if (i > 5):
                                    delta_5 = relative_variations_masses[i_5]
                                    monit["counter"] += 1
                                    if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                        given_mass_dictionary[maximum_mass] += 1
                                        added_here += 1
                                    if (mass_of_new_element > maximum_mass):
                                        
                                        print("break by mass parameter = ", mass_of_new_element)
                                        break
                                for i_6 in range(0,length_relative):
                                    if i <= 6:
                                        break
                                    if (i > 6):
                                        delta_6 = relative_variations_masses[i_6]
                                        monit["counter"] += 1
                                        if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6) == maximum_mass):
                                            given_mass_dictionary[maximum_mass] += 1
                                            added_here += 1
                                        if (mass_of_new_element > maximum_mass):
                                            
                                            print("break by mass parameter = ", mass_of_new_element)
                                            break
                                    for i_7 in range(0,length_relative):
                                        if i <= 7:
                                            break
                                        if (i > 7):
                                            delta_7 = relative_variations_masses[i_7]
                                            monit["counter"] += 1
                                            if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7) == maximum_mass):
                                                given_mass_dictionary[maximum_mass] += 1
                                                added_here += 1
                                            if (mass_of_new_element > maximum_mass):
                                                
                                                print("break by mass parameter = ", mass_of_new_element)
                                                break
                                        for i_8 in range(0,length_relative):
                                            if i <= 8:
                                                break
                                            if (i > 8):
                                                delta_8 = relative_variations_masses[i_8]
                                                monit["counter"] += 1
                                                if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8) == maximum_mass):
                                                    given_mass_dictionary[maximum_mass] += 1
                                                    added_here += 1
                                                if (mass_of_new_element > maximum_mass):
                                                    
                                                    print("break by mass parameter = ", mass_of_new_element)
                                                    break
                                            for i_9 in range(0,length_relative):
                                                if i <= 9:
                                                    break
                                                if (i > 9):
                                                    delta_9 = relative_variations_masses[i_9]
                                                    monit["counter"] += 1
                                                    if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9) == maximum_mass):
                                                        given_mass_dictionary[maximum_mass] += 1
                                                        added_here += 1
                                                    if (mass_of_new_element > maximum_mass):
                                                        
                                                        print("break by mass parameter = ", mass_of_new_element)
                                                        break
                                                for i_10 in range(0,length_relative):
                                                    if i <= 10:
                                                        break
                                                    if (i > 10):
                                                        delta_10 = relative_variations_masses[i_10]
                                                        monit["counter"] += 1
                                                        if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10) == maximum_mass):
                                                            given_mass_dictionary[maximum_mass] += 1
                                                            added_here += 1
                                                        if (mass_of_new_element > maximum_mass):
                                                            
                                                            print("break by mass parameter = ", mass_of_new_element)
                                                            break
                                                    for i_11 in range(0,length_relative):
                                                        if i <= 11:
                                                            break
                                                        if (i > 11):
                                                            delta_11 = relative_variations_masses[i_11]
                                                            monit["counter"] += 1
                                                            if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11) == maximum_mass):
                                                                given_mass_dictionary[maximum_mass] += 1
                                                                added_here += 1
                                                            if (mass_of_new_element > maximum_mass):
                                                                
                                                                print("break by mass parameter = ", mass_of_new_element)
                                                                break
                                                        for i_12 in range(0,length_relative):
                                                            if i <= 12:
                                                                break
                                                            if (i > 12):
                                                                delta_12 = relative_variations_masses[i_12]
                                                                monit["counter"] += 1
                                                                if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12) == maximum_mass):
                                                                    given_mass_dictionary[maximum_mass] += 1
                                                                    added_here += 1
                                                                if (mass_of_new_element > maximum_mass):
                                                                    
                                                                    print("break by mass parameter = ", mass_of_new_element)
                                                                    break
                                                            for i_13 in range(0,length_relative):
                                                                if i <= 13:
                                                                    break
                                                                if (i > 13):
                                                                    delta_13 = relative_variations_masses[i_13]
                                                                    monit["counter"] += 1
                                                                    if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12 + delta_13) == maximum_mass):
                                                                        given_mass_dictionary[maximum_mass] += 1
                                                                        added_here += 1
                                                                    if (mass_of_new_element > maximum_mass):
                                                                        
                                                                        print("break by mass parameter = ", mass_of_new_element)
                                                                        break
                                                                for i_14 in range(0,length_relative):
                                                                    if i <= 14:
                                                                        break
                                                                    if (i > 14):
                                                                        delta_14 = relative_variations_masses[i_14]
                                                                        monit["counter"] += 1
                                                                        if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12 + delta_13 + delta_14) == maximum_mass):
                                                                            given_mass_dictionary[maximum_mass] += 1
                                                                            added_here += 1
                                                                        if (mass_of_new_element > maximum_mass):
                                                                            
                                                                            print("break by mass parameter = ", mass_of_new_element)
                                                                            break
                                                                    for i_15 in range(0,length_relative):
                                                                        if i <= 15:
                                                                            break
                                                                        if (i > 15):
                                                                            delta_15 = relative_variations_masses[i_15]
                                                                            monit["counter"] += 1
                                                                            if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12 + delta_13 + delta_14 + delta_15) == maximum_mass):
                                                                                given_mass_dictionary[maximum_mass] += 1
                                                                                added_here += 1
                                                                            if (mass_of_new_element > maximum_mass):
                                                                                
                                                                                print("break by mass parameter = ", mass_of_new_element)
                                                                                break
                                                                        for i_16 in range(0,length_relative):
                                                                            if i <= 16:
                                                                                break
                                                                            if (i > 16):
                                                                                delta_16 = relative_variations_masses[i_16]
                                                                                monit["counter"] += 1
                                                                                if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12 + delta_13 + delta_14 + delta_15 + delta_16) == maximum_mass):
                                                                                    given_mass_dictionary[maximum_mass] += 1
                                                                                    added_here += 1
                                                                                if (mass_of_new_element > maximum_mass):
                                                                                    
                                                                                    print("break by mass parameter = ", mass_of_new_element)
                                                                                    break
                                                                            for i_17 in range(0,length_relative):
                                                                                if i <= 17:
                                                                                    break
                                                                                if (i > 17):
                                                                                    delta_17 = relative_variations_masses[i_17]
                                                                                    monit["counter"] += 1
                                                                                    if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12 + delta_13 + delta_14 + delta_15 + delta_16 + delta_17) == maximum_mass):
                                                                                        given_mass_dictionary[maximum_mass] += 1
                                                                                        added_here += 1
                                                                                    if (mass_of_new_element > maximum_mass):
                                                                                        
                                                                                        print("break by mass parameter = ", mass_of_new_element)
                                                                                        break
                                                                                for i_18 in range(0,length_relative):
                                                                                    if i <= 18:
                                                                                        break
                                                                                    if (i > 18):
                                                                                        delta_18 = relative_variations_masses[i_18]
                                                                                        monit["counter"] += 1
                                                                                        if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12 + delta_13 + delta_14 + delta_15 + delta_16 + delta_17 + delta_18) == maximum_mass):
                                                                                            given_mass_dictionary[maximum_mass] += 1
                                                                                            added_here += 1
                                                                                        if (mass_of_new_element > maximum_mass):
                                                                                            
                                                                                            print("break by mass parameter = ", mass_of_new_element)
                                                                                            break
                                                                                    for i_19 in range(0,length_relative):
                                                                                        if i <= 19:
                                                                                            break
                                                                                        if (i > 19):
                                                                                            delta_19 = relative_variations_masses[i_19]
                                                                                            monit["counter"] += 1
                                                                                            if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12 + delta_13 + delta_14 + delta_15 + delta_16 + delta_17 + delta_18 + delta_19) == maximum_mass):
                                                                                                given_mass_dictionary[maximum_mass] += 1
                                                                                                added_here += 1
                                                                                            if (mass_of_new_element > maximum_mass):
                                                                                                
                                                                                                print("break by mass parameter = ", mass_of_new_element)
                                                                                                break
                                                                                        for i_20 in range(0,length_relative):
                                                                                            if i <= 20:
                                                                                                break
                                                                                            if (i > 20):
                                                                                                delta_20 = relative_variations_masses[i_20]
                                                                                                monit["counter"] += 1
                                                                                                if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5 + delta_6 + delta_7 + delta_8 + delta_9 + delta_10 + delta_11 + delta_12 + delta_13 + delta_14 + delta_15 + delta_16 + delta_17 + delta_18 + delta_19 + delta_20) == maximum_mass):
                                                                                                    given_mass_dictionary[maximum_mass] += 1
                                                                                                    added_here += 1
                                                                                                if (mass_of_new_element > maximum_mass):
                                                                                                    
                                                                                                    print("break by mass parameter = ", mass_of_new_element)
                                                                                                    break
                                                                                                """
                                                                                            for i_21 in range(0,length_relative):
                                                                                                if i <= 21:
                                                                                                    break
                                                                                                if (i > 21):
                                                                                                    delta_5 = relative_variations_masses[i_21]
                                                                                                    monit["counter"] += 1
                                                                                                    if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                        given_mass_dictionary[maximum_mass] += 1
                                                                                                        added_here += 1
                                                                                                    if (mass_of_new_element > maximum_mass):
                                                                                                        
                                                                                                        print("break by mass parameter = ", mass_of_new_element)
                                                                                                        break
                                                                                                for i_22 in range(0,length_relative):
                                                                                                    if i <= 22:
                                                                                                        break
                                                                                                    if (i > 22):
                                                                                                        delta_5 = relative_variations_masses[i_22]
                                                                                                        monit["counter"] += 1
                                                                                                        if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                            given_mass_dictionary[maximum_mass] += 1
                                                                                                            added_here += 1
                                                                                                        if (mass_of_new_element > maximum_mass):
                                                                                                            
                                                                                                            print("break by mass parameter = ", mass_of_new_element)
                                                                                                            break
                                                                                                    for i_23 in range(0,length_relative):
                                                                                                        if i <= 23:
                                                                                                            break
                                                                                                        if (i > 23):
                                                                                                            delta_5 = relative_variations_masses[i_23]
                                                                                                            monit["counter"] += 1
                                                                                                            if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                                given_mass_dictionary[maximum_mass] += 1
                                                                                                                added_here += 1
                                                                                                            if (mass_of_new_element > maximum_mass):
                                                                                                                
                                                                                                                print("break by mass parameter = ", mass_of_new_element)
                                                                                                                break
                                                                                                        for i_24 in range(0,length_relative):
                                                                                                            if i <= 24:
                                                                                                                break
                                                                                                            if (i > 24):
                                                                                                                delta_5 = relative_variations_masses[i_24]
                                                                                                                monit["counter"] += 1
                                                                                                                if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                                    given_mass_dictionary[maximum_mass] += 1
                                                                                                                    added_here += 1
                                                                                                                if (mass_of_new_element > maximum_mass):
                                                                                                                    
                                                                                                                    print("break by mass parameter = ", mass_of_new_element)
                                                                                                                    break
                                                                                                            for i_25 in range(0,length_relative):
                                                                                                                if i <= 25:
                                                                                                                    break
                                                                                                                if (i > 25):
                                                                                                                    delta_5 = relative_variations_masses[i_25]
                                                                                                                    monit["counter"] += 1
                                                                                                                    if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                                        given_mass_dictionary[maximum_mass] += 1
                                                                                                                        added_here += 1
                                                                                                                    if (mass_of_new_element > maximum_mass):
                                                                                                                        
                                                                                                                        print("break by mass parameter = ", mass_of_new_element)
                                                                                                                        break
                                                                                                                for i_26 in range(0,length_relative):
                                                                                                                    if i <= 26:
                                                                                                                        break
                                                                                                                    if (i > 26):
                                                                                                                        delta_5 = relative_variations_masses[i_26]
                                                                                                                        monit["counter"] += 1
                                                                                                                        if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                                            given_mass_dictionary[maximum_mass] += 1
                                                                                                                            added_here += 1
                                                                                                                        if (mass_of_new_element > maximum_mass):
                                                                                                                            
                                                                                                                            print("break by mass parameter = ", mass_of_new_element)
                                                                                                                            break
                                                                                                                    for i_27 in range(0,length_relative):
                                                                                                                        if i <= 27:
                                                                                                                            break
                                                                                                                        if (i > 27):
                                                                                                                            delta_5 = relative_variations_masses[i_27]
                                                                                                                            monit["counter"] += 1
                                                                                                                            if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                                                given_mass_dictionary[maximum_mass] += 1
                                                                                                                                added_here += 1
                                                                                                                            if (mass_of_new_element > maximum_mass):
                                                                                                                                
                                                                                                                                print("break by mass parameter = ", mass_of_new_element)
                                                                                                                                break
                                                                                                                        for i_28 in range(0,length_relative):
                                                                                                                            if i <= 28:
                                                                                                                                break
                                                                                                                            if (i > 28):
                                                                                                                                delta_5 = relative_variations_masses[i_28]
                                                                                                                                monit["counter"] += 1
                                                                                                                                if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                                                    given_mass_dictionary[maximum_mass] += 1
                                                                                                                                    added_here += 1
                                                                                                                                if (mass_of_new_element > maximum_mass):
                                                                                                                                    
                                                                                                                                    print("break by mass parameter = ", mass_of_new_element)
                                                                                                                                    break
                                                                                                                            for i_29 in range(0,length_relative):
                                                                                                                                if i <= 29:
                                                                                                                                    break
                                                                                                                                if (i > 29):
                                                                                                                                    delta_5 = relative_variations_masses[i_29]
                                                                                                                                    monit["counter"] += 1
                                                                                                                                    if ((mass_of_new_element + delta_1 + delta_2 + delta_3 + delta_4 + delta_5) == maximum_mass):
                                                                                                                                        given_mass_dictionary[maximum_mass] += 1
                                                                                                                                        added_here += 1
                                                                                                                                    if (mass_of_new_element > maximum_mass):
                                                                                                                                        
                                                                                                                                        print("break by mass parameter = ", mass_of_new_element)
                                                                                                                                        break
"""
            print("added_here = ", added_here, "given_mass before = ", given_mass_dictionary[maximum_mass])
            given_mass_dictionary[maximum_mass] += added_here * (i - 1)
            print("given_mass after = ", given_mass_dictionary[maximum_mass])

            return
            
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Submit tasks to the executor
            futures = [executor.submit(adding_other_positions, input_masses) for input_masses in starting_variant_nucleotide_as_a_dict[starting_variant_nucletide]]
            # Collect the results
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
            #print("results = ", results)
        
        print("given_mass_dictionary[maximum_mass] = ", given_mass_dictionary[maximum_mass])
                    
        i += 1

        print("current running time = ", (time.time() - t1))
        
    
    return given_mass_dictionary[maximum_mass]


number_of_variants = counter_of_peptide_variants(1024)

print("number_of_variants = ", number_of_variants)

print("running time = ", (time.time() - t1))