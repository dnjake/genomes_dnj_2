# -*- coding: utf-8 -*-


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np


codes = [('A', 0, 'T', 'U'), ('C', 1, 'G', 'G'), ('G', 2, 'C', 'C'), ('T', 3, 'A', 'A')]

code_dtype = np.dtype([('letter', 'S1'), ('number', np.uint8), ('c_letter', 'S1'), ('c_rna_letter', 'S1')])

codes_a = np.array(codes, code_dtype)
numbers_to_letters = codes_a['letter']
numbers_to_c_letters = codes_a['c_letter']
numbers_to_rna_letters = codes_a['c_rna_letter']

comps_a = np.array([3, 2, 1, 0], dtype=np.uint32)

'''
cpg = 4 + 2
cpg_high = cpg*(2**28)
tg = 6 + 2
tg_high = tg*(2**28)
ca = 5
hi_4 = 15*(2**28)
cg 01100000
96
cg 0110  6
tg 1110  E
ca 0100  4
gg 1010  B
'''


cg_hi = 0X60000000
tg_hi = 0XE0000000
ca_hi = 0X40000000
gg_hi = 0XA0000000
hi_4 =  0XF0000000

def cg_high_mask(nums) :
    vals = np.bitwise_and(nums, hi_4)
    return vals == cg_hi

def tg_high_mask(nums) :
    vals = np.bitwise_and(nums, hi_4)
    return vals == tg_hi

def ca_high_mask(nums) :
    vals = np.bitwise_and(nums, hi_4)
    return vals == ca_hi

def gg_high_mask(nums) :
    vals = np.bitwise_and(nums, hi_4)
    return vals == gg_hi
    
def high_vals(nums) :
    vals = np.bitwise_and(nums, hi_4)
    return vals
    
def or_cg_tg_ca_high_mask(nums) :
    vals = high_vals(nums)
    m = vals == cg_hi
    m = np.logical_or(m, vals == tg_hi)
    m = np.logical_or(m, vals == ca_hi)
    return m

def or_cg_tg_ca_gg_high_mask(nums) :
    vals = high_vals(nums)
    m = vals == cg_hi
    m = np.logical_or(m, vals == tg_hi)
    m = np.logical_or(m, vals == ca_hi)
    m = np.logical_or(m, vals == gg_hi)
    return m


def letters_from_num_16_array(nums) :
    letters_shape = (nums.size, 16)
    letters = np.zeros(letters_shape, dtype='S1')
    for i in range(16) :
        col_numbers = np.bitwise_and(nums, 3)
        col_letters = numbers_to_letters[col_numbers]
        target_col = 15 - i
        letters[:, target_col] = col_letters
        nums = np.right_shift(nums, 2)
    return letters

def num_16_strs(nums) :
    num_16_letters = letters_from_num_16_array(nums)
    out_str = []
    for i in range(nums.size) :
        out_str.append(num_16_letters[i].tostring())
    out_str = np.array(out_str, dtype='S16')
    return out_str

def num_16_short_strs(nums, str_size) :
    num_16_letters = letters_from_num_16_array(nums)
    out_str = []
    for i in range(nums.size) :
        letters = num_16_letters[i]
        out_str.append(letters[:str_size].tostring())
    out_str = np.array(out_str, dtype='S16')
    return out_str

def first_chars_from_num_16_array(nums) :
    top_nums = np.right_shift(nums, 30)
    return numbers_to_letters[top_nums]


def chars_from_num_16(num) :
    out_chars = np.zeros(16, dtype='S1')
    for i in range(16) :
        j = 15 - i
        code = np.bitwise_and(num, 3)
        out_chars[j] = numbers_to_letters[code]
        num = np.right_shift(num, 2)
    return out_chars

def str_from_num_16(num) :
    num_chars = chars_from_num_16(num)
    dna_str = num_chars.tostring()
    return dna_str

def dna_str_from_num_16_array(nums) :
    dna_str = str_from_num_16(nums[-1])
    if nums.size > 1 :
        begin_chars = first_chars_from_num_16_array(nums[:-1])
        dna_str = begin_chars.tostring() + dna_str
    return dna_str
        
def num_16_complement(nums) :
    sel_mask = 3
    out_vals = np.zeros(nums.size, dtype=np.uint32)
    for i in range(16) :
        vals = np.bitwise_and(nums, sel_mask)
        cvals = comps_a[vals]
        out_vals += cvals
        if i < 15 :
            nums = np.right_shift(nums, 2)
            out_vals = np.left_shift(out_vals, 2)
    return out_vals

