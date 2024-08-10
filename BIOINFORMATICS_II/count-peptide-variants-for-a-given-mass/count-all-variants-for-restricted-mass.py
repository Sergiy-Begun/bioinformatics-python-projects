#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 22:13:16 2024

@author: sergiybegun
"""

import time
import copy
import pickle
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
    
    list_of_amino_acid_masses_for_unknown_keys = list(amino_acid_masses_for_unknown.keys())
    
    given_mass_dictionary = {}
    
    length_addons = len(amino_acid_masses_for_unknown.keys())
    
    #thresholds = {}
    
    minimum_addon = min(amino_acid_masses_for_unknown.values())
    
    maximum_addon = max(amino_acid_masses_for_unknown.values())
    
    starting_mass = 2 * maximum_addon
    
    for i in range(starting_mass,(maximum_mass + 1)):
    
        #current_maximum_length_of_peptide = int(float(i) / float(minimum_addon)) + 1
        #print("current_maximum_length_of_peptide = ", current_maximum_length_of_peptide)
        #thresholds[i] = current_maximum_length_of_peptide
        given_mass_dictionary[i] = 0
    
    #current_file_thr = open("thresholds.txt","wb")
    #pickle.dump(thresholds,current_file_thr)
    #current_file_thr.close()
    
    current_file_thr = open("thresholds.txt","rb")
    thr_read_from_file = current_file_thr.read()
    
    thresholds = pickle.loads(thr_read_from_file)
    
    current_file_thr.close()
    
    #print("thresholds = ", thresholds)
    
    #print("thresholds[",maximum_mass,"] = ", thresholds[maximum_mass])
    list_thresholds_keys = list(thresholds.keys())
    
    i = 6
    
    #all_dictionary = copy.deepcopy(amino_acid_masses_for_unknown)
    
    # write for i = 1 only (for other i at the end of cycle)
    #current_file_dict = open(str("all-dictionary-i-" + str(i) + ".txt"),"wb")
    #pickle.dump(all_dictionary,current_file_dict)
    #current_file_dict.close()
    
    current_file_mass = open(str("given_mass_dictionary-i-" + str(i) + ".txt"),"rb")
    given_mass_dictionary_from_file = current_file_mass.read()
    
    given_mass_dictionary = pickle.loads(given_mass_dictionary_from_file)
    current_file_mass.close() 
    
    current_file_dict = open(str("all-dictionary-i-" + str(i) + ".txt"),"rb")
    all_dictionary_from_file = current_file_dict.read()
    
    all_dictionary = pickle.loads(all_dictionary_from_file)
    current_file_dict.close()
    
    print("given_mass_dictionary = ", given_mass_dictionary)
    #print("all_dictionary = ", all_dictionary)
    return
    print("len(all_dictionary.keys()) = ", len(all_dictionary.keys()))

    i = 7
    monit = {"counter": 0}
    while (i < 8):#thresholds[maximum_mass]):
        # formation of peptide dictionary
        print("i = ", i, "len(all_dictionary.keys()) = ", len(all_dictionary.keys()))
        current_level_dictionary = copy.deepcopy(all_dictionary)
        
        total_length = len(all_dictionary.keys())
        
        def adding_new_elements(current_key: str):
            monit["counter"] += 1
            
            if ((monit["counter"] / 100000) == int((monit["counter"] / 100000))):
                print(monit["counter"], "\t", total_length)

            if (len(str(current_key)) < (i - 1)):
                return
        
            for n in range(0,length_addons):
                # print("n = ", n, " of ", length_addons)
                current_new_element = str(current_key) + list_of_amino_acid_masses_for_unknown_keys[n]
                mass_of_new_element = all_dictionary[current_key] + amino_acid_masses_for_unknown[list_of_amino_acid_masses_for_unknown_keys[n]]
                
                current_level_dictionary[current_new_element] = mass_of_new_element
                
                if (mass_of_new_element in given_mass_dictionary.keys()):
                    given_mass_dictionary[mass_of_new_element] += 1
            return
            
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Submit tasks to the executor
            futures = [executor.submit(adding_new_elements, current_key_input) for current_key_input in all_dictionary.keys()]
            # Collect the results
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
            #print("results = ", results)
        
        #print("given_mass_dictionary = ", given_mass_dictionary)
        current_file_mass = open(str("given_mass_dictionary-i-" + str(i) + ".txt"),"wb")
        pickle.dump(given_mass_dictionary,current_file_mass)
        current_file_mass.close()
        
        for thr_i in range(0,len(list_thresholds_keys)):
            #print("list_thresholds_keys[thr_i] = ", list_thresholds_keys[thr_i])
            if (thresholds[list_thresholds_keys[thr_i]] < i):
                graphic_current_file = open(str("graphic.txt"),"a+")
                graphic_current_file.write(str(thr_i) + " = " + str(given_mass_dictionary[list_thresholds_keys[thr_i]]) + "\n")
                graphic_current_file.close()
            
        i += 1
        all_dictionary = copy.deepcopy(current_level_dictionary)
        print("len(all_dictionary.keys() current end of cycle) = ", len(all_dictionary.keys()))
        current_file_dict = open(str("all-dictionary-i-" + str((i - 1)) + ".txt"),"wb")
        pickle.dump(all_dictionary,current_file_dict)
        current_file_dict.close()
        print("current running time = ", (time.time() - t1))
        
    
    return (given_mass_dictionary,all_dictionary)


number_of_variants = counter_of_peptide_variants(1500)

print("len(number_of_variants[0]) = ", len(number_of_variants[0]))

print("len(number_of_variants[1]) = ", len(number_of_variants[1]))

print("running time = ", (time.time() - t1))