"""
useful_funcs.py

This file contains useful functions used throughout the
gPhoton Analysis module. 

Alexander de la Vega (JHU)
Latest version: 5 September 2016
"""

def find_num_harmonic(freq_list, cutoff=0.2):
    count = 0
    n_freqs = len(freq_list)
    if n_freqs != 1:
        for f in freq_list:
            if abs(f - round(f, 0)) <= cutoff:
                if round(f,0) != 0:
                    count += 1
    return count

def find_unique_elements(array):
    found_list = []
    found_list.append(array[0])
    for a in array[1:]:
        if a not in found_list:
            found_list.append(a)
    return found_list

def replace_brackets_and_split(array, char=None):
    if not char:
        return array.replace('[', '').replace(']', '').split()
    else:
        return array.replace('[', '').replace(']', '').split(char)

def check_flag(str_array, flag):
    count = 0
    for j in replace_brackets_and_split(str_array):
        if j == flag:
            count += 1
    if count == 0:
        return False
    else:
        return True