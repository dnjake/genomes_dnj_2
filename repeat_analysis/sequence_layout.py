# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from itertools import izip


class seq_str_cls(object) :
    
    def __init__(self, seq_strs) :
        self.seq_strs = seq_strs
        if type(seq_strs) == str :
            self.seq_count = 1
            self.seq_width = len(self.seq_strs)
        else :
            self.seq_count = len(seq_strs)
            self.seq_width = len(seq_strs[0])
        self.column = None
        
        
class seq_row_cls(object) :
    
    def __init__(self, seq_cols, row_width=80) :
        self.seq_cols = seq_cols
        self.row_strs = None
        self.max_col_count = 1
        self.primary_row = 0
        self.row_width = row_width
    
    def position_cols(self) :
        pos_col = 0
        for seq in self.seq_cols :
            seq.column = pos_col
            pos_col += seq.seq_width + 1
            if seq.seq_count > self.max_col_count :
                self.max_col_count = seq.seq_count

    def write_seq(self, row_str, out_col, out_str) :        
        space_delta = out_col - len(row_str)
        row_str += space_delta*' '
        row_str += out_str
        return row_str
                
    def write_row_strs(self) :
        self.row_strs = []
        for _ in range(self.max_col_count) :
            self.row_strs.append('')        
        for seq in self.seq_cols :
            out_col = seq.column
            if seq.seq_count == 1 :
                out_str = seq.seq_strs
                row = self.primary_row
                row_str = self.row_strs[row]
                self.row_strs[row] = self.write_seq(row_str, out_col, out_str)
            else :
                for i in range(seq.seq_count) :
                    out_str = seq.seq_strs[i]
                    row_str = self.row_strs[i]
                    self.row_strs[i] = self.write_seq(row_str, out_col, out_str)
                    
    def fill_row_ends(self) :
        for i in range(len(self.row_strs)) :
            row_str = self.row_strs[i]
            l = len(row_str)
            if l < self.row_width :
                space_delta = self.row_width - l
                row_str += space_delta*' '
                self.row_strs[i] = row_str
                
    def row_array(self) :
        row_dtype = 'S' + str(self.row_width)
        ra = np.array(self.row_strs, dtype=row_dtype)
        return ra
                
    def do_work(self) :
        self.position_cols()
        self.write_row_strs()
        self.fill_row_ends()
        return self.row_array()
                
                    
class num_16_align_cls(object) :

    def __init__(self, num_counts, letter_group_sizes, pre_space=16, post_space=16) :
        self.num_counts = num_counts
        self.letter_group_sizes = letter_group_sizes
        self.pre_space = pre_space
        self.post_space = post_space
        self.out_string_width = pre_space + post_space + 16 + len(letter_group_sizes) - 1
        self.out_string_dtype = 'S' + str(self.out_string_width)
        self.out_dtype = np.dtype([('num_16', '<u4'), ('count', '<u4'), ('letters', self.out_string_dtype)])
        
    def in_letters_char_view(self) :
        num_letters = self.num_counts['letters'].copy()        
        in_letters = num_letters.view(dtype=('S1',16))
        return in_letters
    
    def build_out_array(self) :
        out_a = np.zeros(self.num_counts.size, dtype=self.out_dtype)
        out_a['num_16'] = self.num_counts['num_16']
        out_a['count'] = self.num_counts['count']
        self.out_num_counts = out_a
        
    def fill_out_array(self) :
        in_letters = self.in_letters_char_view()
        out_shape = (self.num_counts.size, self.out_string_width)
        out_letters = np.zeros(out_shape, dtype='S1')
        out_letters[:] = ' '
        in_col = 0
        out_col = self.pre_space
        for group_size in self.letter_group_sizes :
            in_bound = in_col + group_size
            out_bound = out_col + group_size
            out_letters[:, out_col:out_bound] = in_letters[:, in_col:in_bound]
            in_col = in_bound
            out_col = out_bound + 1
        
        out_strings = out_letters.view(dtype=self.out_string_dtype)
        osr = out_strings.reshape(out_shape[0])
        self.out_num_counts['letters'] = osr
        
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
        
        
        
                