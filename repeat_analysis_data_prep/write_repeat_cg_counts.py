# -*- coding: utf-8 -*-


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import numpy as np
import tables as tb

import genomes_dnj_2.repeat_analysis_data.alu_16_base_data as abd
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d
filters = tb.Filters(complevel=5, complib='zlib')

class data_writer_cls(object) :
    alu_sequence_data_folder = 'grch37_hg19_alu_data'
    cg_offset_dtype = np.dtype([('repeat_index', np.uint32), ('alu_offset', np.uint16)])
    offset_count_dtype = np.dtype([('alu_offset', np.uint16), ('repeat_count', np.uint32), ('cg_count', np.uint32),
                                   ('tg_count', np.uint32), ('ca_count', np.uint32)])
    repeat_index_range_dtype = np.dtype([('repeat_index', np.uint32), ('start', np.uint32), ('bound', np.uint32)])
    repeat_count_dtype = np.dtype([('repeat_index', np.uint32), ('cg_count', np.uint16)])
    count_dist_dtype = np.dtype([('cg_count', np.uint16), ('repeat_count', np.uint32)])
    def __init__(self, data_obj) :
        self.data_obj = data_obj
        
    def build_offset_data(self) :
        d = self.data_obj.num_data
        m = n16d.cg_high_mask(d['num_16'])
        cg_data_size = m.sum()
        cg_array = np.zeros(cg_data_size, dtype=self.cg_offset_dtype)
        cg_array['repeat_index'] = d['repeat_index'][m]
        cg_array['alu_offset'] = d['alu_offset'][m]
        self.cg_repeat_offsets = cg_array

    '''        
    def cg_count_by_offset(self) :
        self.cg_repeat_offsets.sort(order=['alu_offset', 'repeat_index'])
        offsets, counts = np.unique(self.cg_repeat_offsets['alu_offset'], return_counts=True)
        cga = np.zeros(offsets.size, self.offset_count_dtype)
        cga['alu_offset'] = offsets
        cga['cg_count'] = counts
        self.cg_offset_counts = cga
    '''
        
    def cg_count_by_repeat(self) :
        self.cg_repeat_offsets.sort(order=['repeat_index', 'alu_offset'])
        ris = self.data_obj.repeats['index']        
        ri_ranges = np.zeros(ris.size, dtype=self.repeat_index_range_dtype)
        ri_ranges['repeat_index'] = ris
        dris = self.cg_repeat_offsets['repeat_index']
        ri_ranges['start'] = dris.searchsorted(ris)
        ris = ris + 1
        ri_ranges['bound'] = dris.searchsorted(ris)
        out_data = []
        for repeat_index, start, bound in ri_ranges :
            cg_count = bound - start
            out_data.append((repeat_index, cg_count))
        self.repeat_cg_counts = np.array(out_data, dtype=self.repeat_count_dtype)
        
    def cg_repeat_count_distribution(self) :
        cg_counts, repeat_counts = np.unique(self.repeat_cg_counts['cg_count'], return_counts=True)
        da = np.zeros(cg_counts.size, dtype=self.count_dist_dtype)
        da['cg_count'] = cg_counts
        da['repeat_count'] = repeat_counts
        self.cg_count_dist = da
        
    def merge_repeat_data(self) :
        repeats = self.data_obj.repeats
        repeat_dtype = repeats.dtype
        merged_repeat_descr = repeat_dtype.descr
        merged_repeat_descr.append(('cg_count', np.uint16))
        merged_repeats = np.zeros(repeats.size, dtype=np.dtype(merged_repeat_descr))
        repeat_field_names = repeat_dtype.names
        for name in repeat_field_names :
            merged_repeats[name] = repeats[name]
        ris = self.repeat_cg_counts['repeat_index']
        inds_ris = merged_repeats['index'].searchsorted(ris)
        merged_repeats['cg_count'][inds_ris] = self.repeat_cg_counts['cg_count']
        self.merged_repeats = merged_repeats
        
    def do_offset_counts(self) :
        d = self.data_obj.num_data
        d.sort(order='alu_offset')
        offsets, starts, counts = np.unique(d['alu_offset'], return_index=True, return_counts=True)
        out_data = []
        for offset, start, count in zip(offsets, starts, counts) :
            bound = start + count
            od = d[start:bound]
            nums = od['num_16']
            repeat_count = nums.size            
            m = n16d.cg_high_mask(nums)
            cg_count = m.sum()
            m = n16d.tg_high_mask(nums)
            tg_count = m.sum()
            m = n16d.ca_high_mask(nums)
            ca_count = m.sum()
            data = (offset, repeat_count, cg_count, tg_count, ca_count)
            out_data.append(data)
        self.cg_offset_counts = np.array(out_data, dtype=self.offset_count_dtype)
        
        
    def write_output(self) :
        file_path = os.path.join(self.alu_sequence_data_folder, self.file_name)
        h5 = tb.open_file(file_path, 'w', filters=filters)
        index_offset_table = h5.create_table('/', self.repeat_index_cg_offset_table_name,
                                             description=self.cg_repeat_offsets.dtype)
        index_offset_table.append(self.cg_repeat_offsets)
        offset_count_table = h5.create_table('/', self.offset_cg_count_table_name,
                                             description=self.cg_offset_counts.dtype)
        offset_count_table.append(self.cg_offset_counts)
        dist_repeat_cg_count_table = h5.create_table('/', self.dist_repeat_cg_counts_table_name,
                                                     description=self.cg_count_dist.dtype)
        dist_repeat_cg_count_table.append(self.cg_count_dist)
        merged_repeat_table = h5.create_table('/', self.alu_repeat_table_name, 
                                              description=self.merged_repeats.dtype)
        merged_repeat_table.append(self.merged_repeats)
        h5.close()
        
    def do_work(self) :
        self.build_offset_data()        
        self.cg_count_by_repeat()
        self.cg_repeat_count_distribution()
        self.merge_repeat_data()
        self.do_offset_counts()
        self.write_output()

class pos_data_writer_cls(data_writer_cls) :
    file_name = 'pos_strand_repeat_cg_data.h5'    
    alu_repeat_table_name = 'alu_plus_strand_repeats'
    repeat_index_cg_offset_table_name = 'pos_repeat_cg_offsets'
    offset_cg_count_table_name = 'pos_offset_cg_counts'
    dist_repeat_cg_counts_table_name = 'pos_repeat_cg_count_distribution'
    

    def __init__(self) :
        data_obj = abd.pos_alu_data_cls()
        data_writer_cls.__init__(self, data_obj)
        

class neg_data_writer_cls(data_writer_cls) :
    file_name = 'neg_strand_repeat_cg_data.h5'    
    alu_repeat_table_name = 'alu_minus_strand_repeats'
    repeat_index_cg_offset_table_name = 'neg_repeat_cg_offsets'
    offset_cg_count_table_name = 'neg_offset_cg_counts'
    dist_repeat_cg_counts_table_name = 'neg_repeat_cg_count_distribution'
    

    def __init__(self) :
        data_obj = abd.neg_alu_data_cls()
        data_writer_cls.__init__(self, data_obj)

       
dwo = neg_data_writer_cls()

'''
dwo.do_work()
'''
    
    
