# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d

class by_offset_seq_map_cls(object) :
    offset_map_dtype = np.dtype([('alu_offset', np.uint16), ('map_index', np.uint16)])
    repeat_index_map_dtype = np.dtype([('repeat_index', np.uint32), ('map_index', np.uint32)])
    offset_count_dtype = np.dtype([('alu_offset', np.uint16), ('count', np.uint32)])
    repeat_count_dtype = np.dtype([('repeat_index', np.uint32), ('count', np.uint32)])
    def __init__(self, offsets, num_data_obj) :
        self.offset_map = np.zeros(offsets.size, dtype=self.offset_map_dtype)
        self.offset_map['alu_offset'] = offsets
        self.offset_map['map_index'] = np.arange(self.offset_map.size, dtype=np.uint16)
        self.offsets = self.offset_map['alu_offset']
        self.num_data_obj = num_data_obj
        repeat_indexes = self.num_data_obj.repeats['index']
        self.repeat_index_map = np.zeros(repeat_indexes.size, dtype=self.repeat_index_map_dtype)
        self.repeat_index_map['repeat_index'] = repeat_indexes
        self.repeat_index_map['map_index'] = np.arange(repeat_indexes.size, dtype=np.uint32)
        self.repeat_indexes = self.repeat_index_map['repeat_index']
        seq_map_shape = repeat_indexes.size, offsets.size
        self.seq_map = np.zeros(seq_map_shape, dtype=np.uint8)
        
    def build_map(self) :
        for offset, offset_map_index in self.offset_map :
            match_ris = self.find_offset_matches(offset)
            ri_map_indexes = self.repeat_indexes.searchsorted(match_ris)
            self.seq_map[ri_map_indexes, offset_map_index] = 1 

    
    def get_offset_counts(self) :
        cd = np.zeros(self.offset_map.size, dtype=self.offset_count_dtype)
        cd['alu_offset'] = self.offset_map['alu_offset']
        cd['count'] = self.seq_map.sum(axis=0)
        return cd
        
    def get_repeat_counts(self) :
        rd = np.zeros(self.repeat_index_map.size, dtype=self.repeat_count_dtype)
        rd['repeat_index'] = self.repeat_index_map['repeat_index']
        rd['count'] = self.seq_map.sum(axis=1)
        return rd
    
    def get_repeats_indexes_from_offset(self, offset) :
        offset_index = self.offset_map['alu_offset'].searchsorted(offset)
        m = self.seq_map[:, offset_index]
        m = m.astype(np.bool)
        ris = self.repeat_index_map['repeat_index'][m]
        return ris
    
class cg_derived_by_offset_seq_map_cls(by_offset_seq_map_cls) :
    
    def __init__(self, offsets, num_data_obj) :
        by_offset_seq_map_cls.__init__(self, offsets, num_data_obj)    
    
    def find_offset_matches(self, offset) :
        num_data = self.num_data_obj.num_data_from_offset(offset)
        m = n16d.or_cg_tg_ca_high_mask(num_data['num_16'])
        repeat_indexes = num_data['repeat_index'][m]
        return repeat_indexes

class cg_by_offset_seq_map_cls(by_offset_seq_map_cls) :
    
    def __init__(self, offsets, num_data_obj) :
        by_offset_seq_map_cls.__init__(self, offsets, num_data_obj)    
    
    def find_offset_matches(self, offset) :
        num_data = self.num_data_obj.num_data_from_offset(offset)
        m = n16d.cg_high_mask(num_data['num_16'])
        repeat_indexes = num_data['repeat_index'][m]
        return repeat_indexes








