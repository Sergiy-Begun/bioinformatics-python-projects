#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 16:52:49 2024

@author: sergiybegun
"""

import copy
import time

t1 = time.time()


def cyclopeptide_reconstruction_by_convolution_matrix(input_spectrum: list):
    """
    

    Parameters
    ----------
    input_spectrum : list
        DESCRIPTION.

    Returns
    -------
    Spectrum of most frequent elements.

    """
    
    output_spectrum = []
    
    convolution_matrix = []
    
    spectrum_length = len(input_spectrum)
    # print("spectrum_length = ", spectrum_length)
    
    for i in range(0,(spectrum_length - 1)):
        convolution_matrix.append([])
        for j in range((i + 1),spectrum_length):
            convolution_matrix[i].append((input_spectrum[j] - input_spectrum[i]))
    
    # print("convolution_matrix = ", convolution_matrix)
    
    for i in range(0,(spectrum_length - 1)):
        for j in range(0,(spectrum_length - i - 1)):
            # print("i = ", i, "j = ", j)
            current_convolution_matrix_element = convolution_matrix[i][j]
            # print("current_convolution_matrix_element[", i," , ", j, "] = ", current_convolution_matrix_element)
            if (current_convolution_matrix_element > 0):
                output_spectrum.append(current_convolution_matrix_element)
    
    output_spectrum = copy.deepcopy(sorted(output_spectrum,reverse=False))
    
    # print("output_spectrum = ", output_spectrum)
    
    return output_spectrum


read_data_from_file = open("dataset_30246_4.txt", "r")

read_strings_from_file = read_data_from_file.read()

spectrum_list = str(read_strings_from_file).split()

for m in range(0,len(spectrum_list)):
    spectrum_list[m] = int(str(spectrum_list[m]).strip())

spectrum_list = copy.deepcopy(sorted(spectrum_list))

print("spectrum_list = ", spectrum_list)

reconstructed_spectrum = cyclopeptide_reconstruction_by_convolution_matrix(spectrum_list)

output_file = open("output_reconstructed_spectrum.txt", "w")

output_file.write(str(reconstructed_spectrum).replace("\"", "").replace("\'","").replace("\"", "").replace("\'","").replace(",","").replace("[", "").replace("]", ""))

output_file.close()

read_data_from_file.close()

print(time.time() - t1)
