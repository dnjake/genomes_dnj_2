# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')

import genomes_dnj_2.repeat_analysis_data.alu_16_base_data as abd



class count_num_16s_cls(object) :    
    data_dtype = np.dtype([('num_16', np.uint32), ('total_count', np.uint32), ('repeat_count', np.uint32)])
    
    def __init__(self, num_data) :
        self.num_data = num_data
        
    def count_nums(self) :
        self.num_data.sort(order=['num_16', 'repeat_index'])
        nums, starts, counts = np.unique(self.num_data['num_16'], return_index=True, return_counts=True)
        self.num_counts = np.zeros(nums.size, dtype=self.data_dtype)
        self.num_counts['num_16'] = nums
        self.num_counts['total_count'] = counts
        repeat_counts = self.num_counts['repeat_count']
        for i in xrange(nums.size) :
            start = starts[i]
            count = counts[i]
            bound = start + count
            repeat_indexes = self.num_data['repeat_index'][start:bound]
            repeat_indexes = np.unique(repeat_indexes)
            repeat_counts[i] = repeat_indexes.size
            
class merged_num_16_counts_cls(object) :            
    data_dtype = np.dtype([('num_16', np.uint32), ('pos_total_count', np.uint32), ('pos_repeat_count', np.uint32),
                           ('neg_total_count', np.uint32), ('neg_repeat_count', np.uint32)])
    alu_sequence_data_folder = 'grch37_hg19_alu_data'
    count_file_name = 'merged_num_16_counts.h5'
    file_path = os.path.join(alu_sequence_data_folder, count_file_name)
    count_table_name = 'merged_num_16_counts'
    pos_count_table_name = 'pos_num_16_counts'
    neg_count_table_name = 'neg_num_16_counts'
    
    
    def do_counts(self) :
        dobj = abd.pos_alu_data_cls()
        cobj = count_num_16s_cls(dobj.num_data)
        cobj.count_nums()
        self.pos_counts = cobj.num_counts
        dobj = abd.neg_alu_data_cls()
        cobj = count_num_16s_cls(dobj.num_data)
        cobj.count_nums()
        self.neg_counts = cobj.num_counts

    def merge_counts(self) :
        pc = self.pos_counts
        nc = self.neg_counts
        pnums = pc['num_16']
        nnums = nc['num_16']
        nums = np.union1d(pnums, nnums)
        self.merged_counts = np.zeros(nums.size, dtype=self.data_dtype)
        mc = self.merged_counts
        pinds = nums.searchsorted(pnums)
        mc['pos_total_count'][pinds] = pc['total_count']
        mc['pos_repeat_count'][pinds] = pc['repeat_count']
        ninds = nums.searchsorted(nnums)
        mc['neg_total_count'][ninds] = nc['total_count']
        mc['neg_repeat_count'][ninds] = nc['repeat_count']
        
    def write_counts(self) :
        h5 = tb.open_file(self.file_path, 'w', filters=filters)
        count_table = h5.create_table('/', self.count_table_name, description=self.data_dtype)
        count_table.append(self.merged_counts)
        pos_count_table = h5.create_table('/', self.pos_count_table_name, description=self.pos_counts.dtype)
        pos_count_table.append(self.pos_counts)
        neg_count_table = h5.create_table('/', self.neg_count_table_name, description=self.neg_counts.dtype)
        neg_count_table.append(self.neg_counts)
        h5.close()
        
    def do_work(self) :
        self.do_counts()
        self.merge_counts()
        self.write_counts()

'''
mno = merged_num_16_counts_cls()
mno.do_work()        
'''        
        
        