# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d
import genomes_dnj_2.repeat_analysis.parital_sequences_analysis as psa


class match_atom(object) :
    def __init__(self, base_offset, match_delta_offset, match_seqs, match_repeat_indexes) :
        self.base_offset = base_offset,
        self.match_delta_offset = match_delta_offset
        self.match_seqs = match_seqs
        self.match_repeat_indexes = match_repeat_indexes
        
class base_offset_match_cls(object) :
    def __init__(self, base_offset, match_atoms, match_repeat_indexes, contrib_repeat_indexes) :
        self.base_offset = base_offset
        self.match_atoms = match_atoms
        self.match_repeat_indexes = match_repeat_indexes
        self.contrib_repeat_indexes = contrib_repeat_indexes

class match_pattern_cls(object) :
    def __init__(self, delta_offset, match_seqs) :
        self.delta_offset = delta_offset
        self.match_seqs = match_seqs
        


class offset_cluster_match_cls(object) :
    
    def __init__(self, base_offsets, patterns, osq_obj) :
        self.base_offsets = base_offsets
        self.patterns = patterns
        self.osq_obj = osq_obj
        self.base_offset_matches = []
        self.match_repeat_indexes = None
        
    def do_base_offset_match(self, base_offset) :
        match_atoms = []
        for p in self.patterns :
            match_offset = base_offset + p.delta_offset
            p_ris = self.osq_obj.ris_from_offset_seqs(match_offset, p.match_seqs)
            match_atoms.append(match_atom(base_offset, p.delta_offset, p.match_seqs, p_ris))
        ma = match_atoms[0]
        match_repeat_indexes = ma.match_repeat_indexes
        for ma in match_atoms[1:] :
            match_repeat_indexes = np.union1d(match_repeat_indexes, ma.match_repeat_indexes)
        if self.match_repeat_indexes is None :
            self.match_repeat_indexes = match_repeat_indexes
            contrib_repeat_indexes = match_repeat_indexes
        else :
            contrib_repeat_indexes = np.setdiff1d(match_repeat_indexes, self.match_repeat_indexes)
            self.match_repeat_indexes = np.union1d(self.match_repeat_indexes, match_repeat_indexes)
        bmo = base_offset_match_cls(base_offset, self.patterns, match_repeat_indexes, contrib_repeat_indexes)
        self.base_offset_matches.append(bmo)
        
    def do_base_offsets(self) :
        for offset in self.base_offsets :
            self.do_base_offset_match(offset)
            
class cluster_match_num_counts_cls(object) :
    num_count_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('letters', 'S16')])
    
    def __init__(self, base_offset_matches, num_data_obj) :
        self.base_offset_matches = base_offset_matches
        self.num_data_obj = num_data_obj
        
    def build_match_num_counts(self, base_offset_delta=0) :
        match_num_data = []
        for bmo in self.base_offset_matches :
            num_offset = bmo.base_offset + base_offset_delta
            ond = self.num_data_obj.num_data_from_offset_and_repeat_indexes(num_offset, bmo.contrib_repeat_indexes)
            match_num_data.append(ond)
        match_num_data = np.concatenate(match_num_data)
        nums, counts = np.unique(match_num_data['num_16'], return_counts=True)
        mnc = np.zeros(nums.size, dtype=self.num_count_dtype)
        mnc['num_16'] = nums
        mnc['count'] = counts
        mnc['letters'] = n16d.num_16_strs(mnc['num_16'])
        self.num_counts = mnc
        
    def build_sorted_match_num_counts(self, base_offset_delta=0) :
        self.build_match_num_counts(base_offset_delta)
        self.num_counts.sort(order='count')
        self.num_counts = self.num_counts[::-1]
                                          
                                          
class cluster_match_seq_counts_cls(object) :
    seq_count_dtype = np.dtype([('letters', 'S16'), ('count', np.uint32)])
    
    def __init__(self, base_offset_matches, num_data_obj) :
        self.base_offset_matches = base_offset_matches
        self.num_data_obj = num_data_obj
        
    def build_match_seq_counts(self, seq_size, base_offset_delta=0) :
        match_nums = []
        for bmo in self.base_offset_matches :
            num_offset = bmo.base_offset + base_offset_delta
            ond = self.num_data_obj.num_data_from_offset_and_repeat_indexes(num_offset, bmo.contrib_repeat_indexes)
            mn = psa.masked_nums(ond['num_16'], seq_size)
            match_nums.append(mn)
        match_nums = np.concatenate(match_nums)
        nums, counts = np.unique(match_nums, return_counts=True)
        mnc = np.zeros(nums.size, dtype=self.seq_count_dtype)
        mnc['letters'] = n16d.num_16_short_strs(nums, seq_size)
        mnc['count'] = counts
        self.seq_counts = mnc
        
    def build_sorted_match_seq_counts(self, seq_size, base_offset_delta=0) :
        self.build_match_seq_counts(seq_size, base_offset_delta)
        self.seq_counts.sort(order='count')
        self.seq_counts = self.seq_counts[::-1]      
                                          
        
class no_match_num_counts(object) :
    num_count_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('letters', 'S16')])
    def __init__(self, match_ris, num_data_obj) :
        self.match_ris = match_ris
        self.num_data_obj = num_data_obj
        all_ris = num_data_obj.data_all_repeat_indexes()
        self.no_match_ris = np.setdiff1d(all_ris, self.match_ris)
        
    def offset_num_counts(self, offset) :
        no_match_num_data = self.num_data_obj.num_data_from_offset_and_repeat_indexes(offset, self.no_match_ris)
        nums, counts = np.unique(no_match_num_data['num_16'], return_counts=True)
        mnc = np.zeros(nums.size, dtype=self.num_count_dtype)
        mnc['num_16'] = nums
        mnc['count'] = counts
        mnc['letters'] = n16d.num_16_strs(mnc['num_16'])
        return mnc

    def sorted_offset_num_counts(self, offset) :
        num_counts = self.offset_num_counts(offset)
        num_counts.sort(order='count')
        num_counts = num_counts[::-1]
        return num_counts
        
            
        
        
