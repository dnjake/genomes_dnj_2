# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

import os
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')

import genomes_dnj_2.repeat_analysis_data.alu_16_base_data as abd


class offset_num_counts_cls(object) :
    alu_sequence_data_folder = 'grch37_hg19_alu_data'
    offset_num_count_dtype = np.dtype([('num_16', np.uint32), ('alu_offset', np.uint16),
                                       ('offset_count', np.uint32), ('total_count', np.uint32)])
    summary_num_count_dtype = np.dtype([('num_16', np.uint32), ('total_count', np.uint32), ('max_offset', np.uint16),
                                        ('max_offset_count', np.uint32), ('cnt_offsets', np.uint16),
                                        ('cnt_offsets_100_1000', np.uint16), ('cnt_offsets_abv_1000', np.uint16)])
    
    index_dtype = np.dtype([('alu_offset', np.uint32), ('start', np.uint32), ('count', np.uint32)])
    num_index_dtype = np.dtype([('num_16', np.uint32), ('start', np.uint32), ('count', np.uint32)])
    summary_bounds = np.array([100, 1000])
    
    def __init__(self, data_obj, max_offset=300) :
        self.data_obj = data_obj
        self.max_offset = max_offset
        self.num_data = self.data_obj.num_data
        
    def sort_nums(self) :
        self.num_data.sort(order=['alu_offset', 'num_16'])
        offsets, starts, counts = np.unique(self.num_data['alu_offset'], return_index=True, return_counts=True)
        odx = np.zeros(offsets.size, dtype=self.index_dtype)
        odx['alu_offset'] = offsets
        odx['start'] = starts
        odx['count'] = counts
        bound_ind = odx['alu_offset'].searchsorted(self.max_offset + 1)        
        self.offset_index = odx[:bound_ind].copy()
        
    def process_offset(self, offset, num_data) :
        nums, counts = np.unique(num_data['num_16'], return_counts=True)
        oncs = np.zeros(nums.size, dtype=self.offset_num_count_dtype)
        oncs['num_16'] = nums
        oncs['alu_offset'] = offset
        oncs['offset_count'] = counts
        return oncs
    
    def process_offsets(self) :
        out_offset_counts = []
        for offset, start, count in self.offset_index :
            bound = start + count
            offset_num_data = self.num_data[start:bound]
            offset_counts = self.process_offset(offset, offset_num_data)
            out_offset_counts.append(offset_counts)
        self.offset_num_counts = np.concatenate(out_offset_counts)
        
    def sort_offset_num_counts(self) :
        self.offset_num_counts.sort(order=['num_16', 'offset_count'])
        nums, starts, counts = np.unique(self.offset_num_counts['num_16'], return_index=True, return_counts=True)
        ndx = np.zeros(nums.size, dtype=self.num_index_dtype)
        ndx['num_16'] = nums
        ndx['start'] = starts
        ndx['count'] = counts
        self.offset_num_counts_index = ndx
        
    def process_summaries(self) :
        out_summary_data = []
        for num_16, start, count in self.offset_num_counts_index :
            bound = start + count
            nod = self.offset_num_counts[start:bound]
            offset_num, max_offset, max_offset_count, zero_total_count = nod[-1]
            num_total_count = nod['offset_count'].sum()
            nod['total_count'] = num_total_count
            cnt_offsets = nod.size
            bound_inds = nod['offset_count'].searchsorted(self.summary_bounds)
            cnt_100_1000 = bound_inds[1] - bound_inds[0]
            cnt_abv_1000 = cnt_offsets - bound_inds[1]
            data = (num_16, num_total_count, max_offset, max_offset_count, cnt_offsets, cnt_100_1000, cnt_abv_1000)
            out_summary_data.append(data)
        self.num_offset_summary_data = np.array(out_summary_data, dtype=self.summary_num_count_dtype)
        
    def write_data(self) :
        file_path = os.path.join(self.alu_sequence_data_folder, self.offset_num_counts_file_name)
        h5 = tb.open_file(file_path, 'w', filters=filters)
        offset_num_counts_table = h5.create_table('/', self.offset_num_counts_table_name,
                                                  description=self.offset_num_counts.dtype)
        offset_num_counts_table.append(self.offset_num_counts)
        summary_num_counts_table = h5.create_table('/', self.summary_num_counts_table_name,
                                                   description=self.num_offset_summary_data.dtype)
        summary_num_counts_table.append(self.num_offset_summary_data)
        h5.close()
        
    def do_work(self) :
        self.sort_nums()
        self.process_offsets()
        self.sort_offset_num_counts()
        self.process_summaries()
        self.write_data()
                

class pos_offset_num_counts_cls(offset_num_counts_cls) :
    offset_num_counts_file_name = 'pos_alu_offset_num_counts.h5'
    offset_num_counts_table_name = 'pos_alu_offset_num_counts'    
    summary_num_counts_table_name = 'pos_alu_num_offset_count_summaries'    
    
    def __init__(self) :
        data_obj = abd.pos_alu_data_cls()
        offset_num_counts_cls.__init__(self, data_obj)
        
        
class neg_offset_num_counts_cls(offset_num_counts_cls) :
    offset_num_counts_file_name = 'neg_alu_offset_num_counts.h5'
    offset_num_counts_table_name = 'neg_alu_offset_num_counts'    
    summary_num_counts_table_name = 'neg_alu_num_offset_count_summaries'    
    
    def __init__(self) :
        data_obj = abd.neg_alu_data_cls()
        offset_num_counts_cls.__init__(self, data_obj)
        
'''        
onco = neg_offset_num_counts_cls()
onco.do_work()        
'''
        
    
        
