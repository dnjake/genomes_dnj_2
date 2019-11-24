# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d
import genomes_dnj_2.repeat_analysis.parital_sequences_analysis as psa


class sequence_value_match_cls(object) :
    match_data_dtype = np.dtype([('repeat_index', np.uint32), ('match_value', np.uint32), ('match_offset', np.uint16),
                                 ('assoc_offset', np.uint16), ('assoc_value', np.uint32), ('assoc_valid', np.bool)])
    offset_match_dtype = np.dtype([('offset', np.uint16), ('count', np.uint32)])
    value_match_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('letters', 'S16')])
    subseq_match_dtype = np.dtype([('seq', 'S16'), ('count', np.uint32)])
    
    def __init__(self, num_data, subseq_match=False, mask_size=0) :
        self.num_data = num_data
        self.subseq_match = subseq_match
        self.mask_size = mask_size
        if self.subseq_match :
            self.num_data['num_16'] = psa.masked_nums(self.num_data['num_16'], self.mask_size)
        self.num_data.sort(order=['num_16', 'repeat_index'])
        
        
    def match_num_16_value(self, match_value, test_repeat_indexes) :
        value_bound = match_value + 1
        inds = self.num_data['num_16'].searchsorted([match_value, value_bound])
        match_data = self.num_data[inds[0]:inds[1]]
        if match_data.size == 0 :
            return None
        inds = match_data['repeat_index'].searchsorted(test_repeat_indexes)
        m = inds < match_data.size
        data_inds = inds[m]
        test_inds = np.where(m)[0] 
        m = match_data['repeat_index'][data_inds] == test_repeat_indexes[test_inds]
        data_size = m.sum()
        if data_size == 0 :
            return None
        data_inds = data_inds[m]
        test_inds = test_inds[m]
        data = np.zeros(data_size, self.match_data_dtype)
        data['repeat_index'] = test_repeat_indexes[test_inds]
        data['match_value'] = match_data['num_16'][data_inds]
        data['match_offset'] = match_data['alu_offset'][data_inds]
        return data
    
    def match_values(self, match_values, test_repeat_indexes) :
        if type(match_values) == int :
            match_values = [match_values]
        matches = []
        for match_value in match_values :
            match_data = self.match_num_16_value(match_value, test_repeat_indexes)
            if match_data is None :
                continue
            matches.append(match_data)
            matched_repeat_indexes = match_data['repeat_index']
            test_repeat_indexes = np.setdiff1d(test_repeat_indexes, matched_repeat_indexes)
        if len(matches) == 0 :
            return None
        match_data = np.concatenate(matches)
        return match_data
    
    def match_subseqs(self, subseqs, test_repeat_indexes) :
        if type(subseqs) == str :
            subseqs = [subseqs]
        match_values = []
        for subseq in subseqs :
            subseq_val = psa.num_16_from_substr(subseq)
            match_values.append(subseq_val)
        return self.match_values(match_values, test_repeat_indexes)
    
    def matched_offsets(self, match_data) :
        offsets, counts = np.unique(match_data['match_offset'], return_counts=True)
        dt = np.zeros(offsets.size, dtype=self.offset_match_dtype)
        dt['offset'] = offsets
        dt['count'] = counts
        return dt
            
    def matched_num_16s(self, match_data) :
        values, counts = np.unique(match_data['match_value'], return_counts=True)
        dt = np.zeros(values.size, dtype = self.value_match_dtype)
        dt['num_16'] = values
        dt['count'] = counts
        dt['letters'] = n16d.num_16_strs(dt['num_16'])
        dt.sort(order='count')
        dt = dt[::-1]
        return dt
        
    def matched_subseqs(self, match_data) :
        values, counts = np.unique(match_data['match_value'], return_counts=True)
        dt = np.zeros(values.size, dtype=self.subseq_match_dtype)
        dt['seq'] = n16d.num_16_short_strs(values, self.mask_size)
        dt['count'] = counts
        dt.sort(order='count')
        dt = dt[::-1]
        return dt
        
class associated_values_cls(object) :
    offset_dtype = np.dtype([('offset', np.uint16), ('start', np.uint32), ('count', np.uint32)])
    value_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('letters', 'S16')])
    subseq_dtype = np.dtype([('seq', 'S16'), ('count', np.uint32)])
    
    def __init__(self, assoc_offset, match_data, num_data) :
        self.num_data = num_data
        self.match_data = match_data
        self.match_data['assoc_offset'] = self.match_data['match_offset'] + assoc_offset
        self.num_data.sort(order=['alu_offset', 'repeat_index'])
        self.match_data.sort(order=['assoc_offset', 'repeat_index'])
        offsets, starts, counts = np.unique(self.match_data['assoc_offset'], return_index=True, return_counts=True)
        omc = np.zeros(offsets.size, dtype=self.offset_dtype)
        omc['offset'] = offsets
        omc['start'] = starts
        omc['count'] = counts
        self.offset_match_counts = omc
        self.fill_in_assoc_values()

    def do_offset_match(self, offset, offset_match_data) :
        bound = offset + 1
        inds = self.num_data['alu_offset'].searchsorted([offset, bound])
        offset_num_data = self.num_data[inds[0]:inds[1]]
        match_repeat_indexes = offset_match_data['repeat_index']
        inds = offset_num_data['repeat_index'].searchsorted(match_repeat_indexes)
        m_in_num_data = inds < offset_num_data.size
        inds_num_data = inds[m_in_num_data]
        inds_matched_repeats = np.where(m_in_num_data)[0]
        m_matched_repeats = offset_num_data['repeat_index'][inds_num_data] == match_repeat_indexes[inds_matched_repeats]
        inds_matched_repeats = inds_matched_repeats[m_matched_repeats]
        inds_num_data = inds_num_data[m_matched_repeats]
        offset_match_data['assoc_value'][inds_matched_repeats] = offset_num_data['num_16'][inds_num_data]
        offset_match_data['assoc_valid'][inds_matched_repeats] = True
        
    def fill_in_assoc_values(self) :
        for offset, start, count, in self.offset_match_counts :
            bound = start + count
            offset_match_data = self.match_data[start:bound]
            self.do_offset_match(offset, offset_match_data)
            
    def valid_assoc_values(self) :
        return self.match_data['assoc_value'][self.match_data['assoc_valid']]
    
    def associated_num_16s(self) :
        values, counts = np.unique(self.valid_assoc_values(), return_counts=True)
        dt = np.zeros(values.size, dtype=self.value_dtype)
        dt['num_16'] = values
        dt['count'] = counts
        dt['letters'] = n16d.num_16_strs(dt['num_16'])
        dt.sort(order='count')
        dt = dt[::-1]
        return dt
    
    def associated_subseqs(self, seq_size) :
        values, counts = np.unique(self.valid_assoc_values(), return_counts=True)
        dt = np.zeros(values.size, dtype=self.subseq_dtype)
        dt['seq'] = n16d.num_16_short_strs(values, seq_size)
        dt['count'] = counts
        dt.sort(order='count')
        dt = dt[::-1]
        return dt
        
