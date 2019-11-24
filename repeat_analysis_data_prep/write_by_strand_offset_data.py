# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')

import genomes_dnj_2.repeat_analysis_data.alu_16_base_data as abd
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d


class offset_anal_cls(object) :
    #count_dtype = np.dtype([('offset', np.uint32), ('start', np.uint32), ('count', np.uint32)])
    num_count_dtype = np.dtype([('num_16', np.uint32), ('start', np.uint32), ('count', np.uint32)])
    offset_dtype = np.dtype([('alu_offset', np.uint16), ('max_num_16', np.uint32),                              
                     ('max_num_offset_count', np.uint32), ('max_num_total_count', np.uint32),('max_num_repeat_count', np.uint32),
                     ('unique_offset_nums', np.uint32), ('total_offset_count', np.uint32), ('cpg_count', np.uint32),
                     ('nuns_abv_10000', np.uint32), ('count_abv_10000', np.uint32), ('nums_1000_10000', np.uint32),
                     ('count_1000_10000', np.uint32), ('nums_100_1000', np.uint32), ('count_100_1000', np.uint32),
                     ('nums_10_100', np.uint32), ('count_10_100', np.uint32), ('nums_less_10', np.uint32),
                     ('count_less_10', np.uint32)])
    num_count_bins = np.array([10, 100, 1000, 10000], dtype=np.uint32)
    
    def __init__(self, num_data, num_counts) :
        self.num_data = num_data
        self.num_counts = num_counts
        
    def process_offset_num_data(self, alu_offset, offset_num_data) :
        offset_nums, starts, counts = np.unique(offset_num_data['num_16'], return_index=True, return_counts=True)
        offset_num_counts = np.array(zip(offset_nums, starts, counts), dtype=self.num_count_dtype)
        offset_num_counts.sort(order=['count', 'num_16'])
        unique_offset_nums_count = offset_num_counts.size
        total_offset_count = offset_num_counts['count'].sum()
        cpg_mask = n16d.cpg_high_mask(offset_nums)
        cpg_counts = counts[cpg_mask]
        cpg_count = cpg_counts.sum()
        bin_inds = offset_num_counts['count'].searchsorted(self.num_count_bins)
        max_num_16, max_start, max_num_offset_count = offset_num_counts[-1]
        max_num_total_count = 0
        max_num_repeat_count = 0
        offset_data = [alu_offset, max_num_16, max_num_offset_count, max_num_total_count, max_num_repeat_count,
                       unique_offset_nums_count, total_offset_count, cpg_count]
        low_ind = bin_inds[-1]
        bin_num_counts = offset_num_counts[low_ind:]
        bin_nums = bin_num_counts.size
        bin_count = bin_num_counts['count'].sum()
        offset_data.extend([bin_nums, bin_count])
        for i in range(1,4) :
            high_ind = bin_inds[-i]
            low_ind = bin_inds[-i-1]
            bin_num_counts = offset_num_counts[low_ind:high_ind]
            bin_nums = bin_num_counts.size
            bin_count = bin_num_counts['count'].sum()
            offset_data.extend([bin_nums, bin_count])
        high_ind = bin_inds[0]
        bin_num_counts = offset_num_counts[:high_ind]
        bin_nums = bin_num_counts.size
        bin_count = bin_num_counts['count'].sum()
        offset_data.extend([bin_nums, bin_count])
        self.repeat_offset_data.append(tuple(offset_data))

    def process_offsets(self) :
        self.num_data.sort(order=['alu_offset', 'num_16'])
        self.repeat_offset_data = []
        offsets, starts, counts = np.unique(self.num_data['alu_offset'], return_index=True, return_counts=True)
        for offset, start, count in zip(offsets, starts, counts) :
            print('processing offset', offset)
            bound = start + count
            offset_num_data = self.num_data[start:bound]
            self.process_offset_num_data(offset, offset_num_data)
        self.repeat_offset_data = np.array(self.repeat_offset_data, dtype=self.offset_dtype)
        nums = self.repeat_offset_data['max_num_16']
        inds_nums = self.num_counts['num_16'].searchsorted(nums)
        self.repeat_offset_data['max_num_total_count'] = self.num_counts['total_count'][inds_nums]
        self.repeat_offset_data['max_num_repeat_count'] = self.num_counts['repeat_count'][inds_nums]


    def do_work(self) :
        self.process_offsets()


class offset_data_writer_cls(object) :
    data_folder = 'grch37_hg19_alu_data'
    file_name = 'alu_strand_offset_data_stats.h5'
    file_path = os.path.join(data_folder, file_name)
    pos_table_name = 'pos_alu_offset_data_stats'
    neg_table_name = 'neg_alu_offset_data_stats'

    def write_pos_data(self, h5) :
        num_data_obj = abd.pos_alu_data_cls()
        num_data = num_data_obj.num_data
        num_counts_obj = abd.pos_alu_num_counts()
        num_counts = num_counts_obj.num_counts
        oao = offset_anal_cls(num_data, num_counts)
        oao.do_work()
        data = oao.repeat_offset_data
        pos_table = h5.create_table('/', self.pos_table_name, description=data.dtype)
        pos_table.append(data)
        pos_table.close()

    def write_neg_data(self, h5) :
        num_data_obj = abd.neg_alu_data_cls()
        num_data = num_data_obj.num_data
        num_counts_obj = abd.neg_alu_num_counts()
        num_counts = num_counts_obj.num_counts
        oao = offset_anal_cls(num_data, num_counts)
        oao.do_work()
        data = oao.repeat_offset_data
        neg_table = h5.create_table('/', self.neg_table_name, description=data.dtype)
        neg_table.append(data)
        neg_table.close()
        
    def do_work(self) :
        h5 = tb.open_file(self.file_path, 'w', filters=filters)
        self.write_pos_data(h5)
        self.write_neg_data(h5)
        h5.close()

'''        
odwo = offset_data_writer_cls()
odwo.do_work()
'''

