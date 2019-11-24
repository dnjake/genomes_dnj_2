# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from . import offset_analysis as oa

class region_analysis_cls(object) :
    repeat_data_dtype = np.dtype([('repeat_index', np.uint32), ('offset', np.uint16),
                           ('num_16', np.uint32), ('num_count', np.uint32)])
    offset_count_dtype = np.dtype([('offset', np.uint16), ('num_16', np.uint32), ('num_count', np.uint32),
                                   ('repeat_count', np.uint32)])
    
    def __init__(self, offset_low, offset_high, data_obj, min_offset_count=100) :
        self.offset_low = offset_low
        self.offset_high = offset_high
        self.data_obj = data_obj
        self.min_offset_count = min_offset_count
        data_repeat_indexes = self.data_obj.repeats['index']
        repeat_data = np.zeros(data_repeat_indexes.size, dtype=self.repeat_data_dtype)
        repeat_data['repeat_index'] = data_repeat_indexes
        self.repeat_data = repeat_data
        self.repeat_indexes = self.repeat_data['repeat_index']
        self.observed_repeat_indexes = None
        
    def process_offset(self, offset) :
        onao = oa.offset_nums_anal_cls(offset, self.data_obj)
        offset_repeat_indexes = onao.num_data['repeat_index']
        if self.observed_repeat_indexes is None :
            self.observed_repeat_indexes = offset_repeat_indexes.copy()
        else :
            self.observed_repeat_indexes = np.union1d(self.observed_repeat_indexes, offset_repeat_indexes)
        offset_num_counts = onao.gen_count_num_data(self.min_offset_count)
        for num_16, count, count_num_data in offset_num_counts :
            inds_repeats = self.repeat_indexes.searchsorted(count_num_data['repeat_index'])
            m = self.repeat_data['num_count'][inds_repeats] < count
            inds_repeats = inds_repeats[m]
            self.repeat_data['offset'][inds_repeats] = offset 
            self.repeat_data['num_16'][inds_repeats] = num_16
            self.repeat_data['num_count'][inds_repeats] = count
            
    def process_offsets(self) :
        for offset in range(self.offset_low, self.offset_high) :
            self.process_offset(offset)
            
    def sumarize_data(self) :
        m = self.repeat_data['offset'] > 0
        self.repeat_data = self.repeat_data[m]
        self.repeat_data.sort(order=['num_16', 'offset'])
        repeat_nums, repeat_starts, repeat_counts = np.unique(self.repeat_data['num_16'], 
                                                              return_index=True, return_counts=True)
        out_data = []
        for repeat_num, repeat_num_start, repeat_num_count in zip(repeat_nums, repeat_starts, repeat_counts) :
            repeat_num_bound = repeat_num_start + repeat_num_count
            repeat_num_data = self.repeat_data[repeat_num_start:repeat_num_bound]
            offsets, offset_starts, offset_counts = np.unique(repeat_num_data['offset'], return_index=True,
                                                              return_counts=True)
            for offset, offset_start, offset_count in zip(offsets, offset_starts, offset_counts) :
                ri, offset, offset_num_16, offset_num_count = repeat_num_data[offset_start]
                out_data.append((offset, offset_num_16, offset_num_count, offset_count))
        self.region_offset_num_counts = np.array(out_data, dtype=self.offset_count_dtype)
        
    def sorted_region_data(self) :
        self.process_offsets()
        self.sumarize_data()
        self.region_offset_num_counts.sort(order=['repeat_count', 'offset'])
        self.region_offset_num_counts = self.region_offset_num_counts[::-1]
        return self.region_offset_num_counts
        
            
            



