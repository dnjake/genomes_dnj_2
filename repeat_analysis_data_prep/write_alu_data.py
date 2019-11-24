# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb

import genomes_dnj_2.repeat_anal_data.all_repeat_data as ard

filters = tb.Filters(complevel=5, complib='zlib')


class chrom_16_base_rdr_cls(object) :
    sequence_numbers_folder = 'grch37_hg19_chrom_dna_8_base_codes'
    sequence_numbers_file_name_end = '_8_base_sequence_numbers_by_pos.h5'
    sequence_numbers_carray_name_end = '_8_base_sequence_numbers_by_pos'
    valid_carray_name_end = '_valid_pos_mask'
    code_offset = 8
    
    def __init__(self, chrom) :
        self.chrom = chrom
        self.chrom_start_str = 'chr_' + str(self.chrom)
        self.read_chrom_input()
        self.create_seq_views()

    def read_chrom_input(self) :
        file_name = self.chrom_start_str + self.sequence_numbers_file_name_end
        file_path = os.path.join(self.sequence_numbers_folder, file_name)
        in_seq_nums_array_name = self.chrom_start_str + self.sequence_numbers_carray_name_end 
        in_valid_mask_name = self.chrom_start_str + self.valid_carray_name_end
        h5 = tb.open_file(file_path, 'r')
        seq_nums_array = getattr(h5.root, in_seq_nums_array_name)
        self.seq_nums = seq_nums_array[:]
        valid_mask_array = getattr(h5.root, in_valid_mask_name)
        self.valid_mask = valid_mask_array[:]
        h5.close()

    def create_seq_views(self) :
        offset = self.code_offset
        self.first_seqs = self.seq_nums[:-offset]
        self.second_seqs = self.seq_nums[offset:]
        first_valid = self.valid_mask[:-offset]
        second_valid = self.valid_mask[offset:]
        self.valid_16_bases = np.logical_and(first_valid, second_valid)
        
        
class alu_dna_cls(object) :
    alu_sequence_data_folder = 'grch37_hg19_alu_data'
    alu_seq_dtype = np.dtype([('chrom', np.uint16), ('alu_offset', np.uint16), ('pos', np.uint32),
                              ('num_16', np.uint32), ('repeat_index', np.uint32)])
    unique_num_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32)])

    def __init__(self, strand='+') :
        self.strand = strand
        self.alu_data_obj = ard.alu_repeat_cls(self.strand)
        
    
    def do_all_chroms(self) :
        all_alu_data = []
        for chrom in range(1,23) :
            print('processing chrom', chrom)
            chrom_alu_data = self.do_chrom_extract(chrom)
            all_alu_data.append(chrom_alu_data)
        self.all_alu_data = np.concatenate(all_alu_data)

    
    def write_data(self) :
        folder_name = self.alu_sequence_data_folder
        file_name = self.alu_sequence_data_file_name
        file_path = os.path.join(folder_name, file_name)
        data_table_name = self.alu_num_16_table_name
        alu_repeat_table_name = self.alu_repeat_table_name
        self.all_alu_data.sort(order=['repeat_index', 'alu_offset'])
        h5 = tb.open_file(file_path, 'w', filters=filters)
        data_table = h5.create_table('/', data_table_name, description=self.all_alu_data.dtype)
        data_table.append(self.all_alu_data)
        alu_repeat_data = self.alu_data_obj.repeat_data
        alu_repeat_table = h5.create_table('/', alu_repeat_table_name, description=alu_repeat_data.dtype)
        alu_repeat_table.append(alu_repeat_data)
        h5.close()


    def do_write_work(self) :
        self.do_all_chroms()
        self.write_data()
        print('done')
        
        
        
        


class pos_alu_dna_cls(alu_dna_cls) :
    alu_sequence_data_file_name = 'pos_alu_dna_num_16.h5'
    alu_num_16_table_name = 'pos_alu_dna_num_16_data'
    alu_repeat_table_name = 'alu_plus_strand_repeats'
    
    def __init__(self) :
        alu_dna_cls.__init__(self, strand='+')
        
    def strand_data(self, chrom, alu_first_seqs, alu_second_seqs, alu_index,
                    alu_start_pos, alu_bound_pos, alu_start_offset) :
        data_size = alu_first_seqs.size - 15
        alu_data = np.zeros(data_size, self.alu_seq_dtype)
        alu_data['chrom'] = chrom
        alu_data['alu_offset'] = np.arange(data_size, dtype=np.uint16)
        alu_data['alu_offset'] += alu_start_offset
        bound_pos = alu_bound_pos - 15
        alu_data['pos'] = np.arange(alu_start_pos, bound_pos, dtype=np.uint32)
        alu_data['repeat_index'] = alu_index
        alu_data['num_16'] = alu_first_seqs[:data_size]
        alu_data['num_16'] *= 2**16
        alu_data['num_16'] += alu_second_seqs[:data_size]
        return alu_data
        
    def do_chrom_extract(self, chrom) :
        chrom_view_obj = chrom_16_base_rdr_cls(chrom)
        chrom_alu_repeats = self.alu_data_obj.chrom_alu_data(chrom)
        alu_starts = chrom_alu_repeats['start_pos']
        alu_bounds = chrom_alu_repeats['end_pos']
        alu_indexes = chrom_alu_repeats['index']
        alu_offsets = chrom_alu_repeats['repeat_start']
        chrom_alu_data = []
        for i in xrange(chrom_alu_repeats.size) :
            alu_index = alu_indexes[i]
            alu_start_pos = alu_starts[i] + 1
            alu_bound_pos = alu_bounds[i] + 1
            if alu_bound_pos - alu_start_pos < 16 :
                continue
            alu_start_offset = alu_offsets[i]
            valid_mask = chrom_view_obj.valid_16_bases[alu_start_pos:alu_bound_pos]
            if not valid_mask.all() :
                continue
            chrom = chrom_view_obj.chrom
            alu_first_seqs = chrom_view_obj.first_seqs[alu_start_pos:alu_bound_pos]
            alu_second_seqs = chrom_view_obj.second_seqs[alu_start_pos:alu_bound_pos]
            alu_data = self.strand_data(chrom, alu_first_seqs, alu_second_seqs, alu_index, alu_start_pos,
                                           alu_bound_pos, alu_start_offset)
            chrom_alu_data.append(alu_data)
        chrom_alu_data = np.concatenate(chrom_alu_data)
        return chrom_alu_data
        
class neg_alu_dna_cls(alu_dna_cls) :
    alu_sequence_data_file_name = 'neg_alu_dna_num_16.h5'
    alu_num_16_table_name = 'neg_alu_dna_num_16_data'
    alu_repeat_table_name = 'alu_minus_strand_repeats'
    base_comps = np.array([3, 2, 1, 0], dtype=np.uint32)
    
    def __init__(self) :
        alu_dna_cls.__init__(self, strand='-')
        
    def num_complements(self, nums) :
        sel_mask = 3
        out_vals = np.zeros(nums.size, dtype=np.uint32)
        for i in range(16) :
            vals = np.bitwise_and(nums, sel_mask)
            cvals = self.base_comps[vals]
            out_vals += cvals
            if i < 15 :
                nums = np.right_shift(nums, 2)
                out_vals = np.left_shift(out_vals, 2)
        return out_vals


    def strand_data(self, chrom, alu_first_seqs, alu_second_seqs, alu_index,
                    alu_start_pos, alu_bound_pos, alu_start_offset) :
        data_size = alu_first_seqs.size - 15
        alu_data = np.zeros(data_size, self.alu_seq_dtype)
        alu_data['chrom'] = chrom
        alu_data['alu_offset'] = np.arange(data_size, dtype=np.uint16)
        alu_data['alu_offset'] = data_size - alu_data['alu_offset'] -1
        
        alu_data['alu_offset'] += alu_start_offset
        start_pos = alu_start_pos + 15
        alu_data['pos'] = np.arange(start_pos, alu_bound_pos, dtype=np.uint32)
        alu_data['repeat_index'] = alu_index
        nums = np.zeros(data_size, np.uint32)
        nums[:] = alu_first_seqs[:data_size]
        nums *= 2**16
        nums += alu_second_seqs[:data_size]
        alu_data['num_16'] = self.num_complements(nums)
        return alu_data

    def do_chrom_extract(self, chrom) :
        chrom_view_obj = chrom_16_base_rdr_cls(chrom)
        chrom_alu_repeats = self.alu_data_obj.chrom_alu_data(chrom)
        alu_starts = chrom_alu_repeats['start_pos'] + 1
        alu_bounds = chrom_alu_repeats['end_pos'] + 1
        alu_indexes = chrom_alu_repeats['index']
        alu_offsets = chrom_alu_repeats['repeat_left']
        chrom_alu_data = []
        for i in xrange(chrom_alu_repeats.size) :
            alu_index = alu_indexes[i]
            alu_start_pos = alu_starts[i]
            alu_bound_pos = alu_bounds[i]
            if alu_bound_pos - alu_start_pos < 16 :
                continue
            alu_start_offset = alu_offsets[i]
            valid_mask = chrom_view_obj.valid_16_bases[alu_start_pos:alu_bound_pos]
            if not valid_mask.all() :
                continue
            chrom = chrom_view_obj.chrom
            alu_first_seqs = chrom_view_obj.first_seqs[alu_start_pos:alu_bound_pos]
            alu_second_seqs = chrom_view_obj.second_seqs[alu_start_pos:alu_bound_pos]
            alu_data = self.strand_data(chrom, alu_first_seqs, alu_second_seqs, alu_index, alu_start_pos,
                                           alu_bound_pos, alu_start_offset)
            chrom_alu_data.append(alu_data)
        chrom_alu_data = np.concatenate(chrom_alu_data)
        return chrom_alu_data

'''
pos = pos_alu_dna_cls()
pos.do_write_work()
'''

