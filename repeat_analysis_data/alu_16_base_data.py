# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)



class repeat_num_data_cls(object) :
    '''
    This class holds the num data for some set of repeats
    It is assumed to be inititalized with that num data sorted
    by repeat and offset and the repeats that contain the data
    '''
    offset_count_dtype = np.dtype([('alu_offset', np.uint32), ('count', np.uint32)])
    num_count_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32)])
    num_count_letters_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('letters', 'S16')])

    def __init__(self, num_data, repeats) :
        self.num_data = num_data
        self.repeats = repeats
        self.repeat_count = self.repeats.size
    
    def num_data_from_repeat_indexes(self, repeat_indexes) :
        bounds = repeat_indexes + 1
        inds_starts = self.num_data['repeat_index'].searchsorted(repeat_indexes)
        inds_bounds = self.num_data['repeat_index'].searchsorted(bounds)
        out_data = []
        for ind_start, ind_bound in zip(inds_starts, inds_bounds) :
            out_data.append(self.num_data[ind_start:ind_bound])
        out_data = np.concatenate(out_data)
        return out_data            
            
    def data_all_repeat_indexes(self) :
        return self.repeats['index']

    def repeat_indexes_from_num(self, num, offset=None) :
        m = self.num_data['num_16'] == num
        if offset is not None :
            mo = self.num_data['alu_offset'] == offset
            m = np.logical_and(m, mo)
        repeat_indexes = self.num_data['repeat_index'][m]
        if offset is None :
            repeat_indexes = np.unique(repeat_indexes)
        return repeat_indexes

    def repeat_indexes_from_num_and_offset_range(self, num, offset_start, offset_bound) :
        m = self.num_data['num_16'] == num
        m = np.logical_and(m, self.num_data['alu_offset'] >= offset_start)
        m = np.logical_and(m, self.num_data['alu_offset'] < offset_bound)
        repeat_indexes = self.num_data['repeat_index'][m]
        repeat_indexes = np.unique(repeat_indexes)
        return repeat_indexes


    def repeat_indexes_from_offset_nums(self, offset_nums) :
        ris = None
        for offset, num in offset_nums :
            oris = self.repeat_indexes_from_num(num, offset=offset)
            if ris is None :
                ris = oris
            else :
                ris = np.intersect1d(ris, oris)
        return ris

    def num_data_from_offset(self, offset) :
        m = self.num_data['alu_offset'] == offset
        return self.num_data[m]

    def num_data_from_offset_and_repeat_indexes(self, offset, repeat_indexes) :
        offset_num_data = self.num_data_from_offset(offset)
        inds_data = offset_num_data['repeat_index'].searchsorted(repeat_indexes)
        m = inds_data < offset_num_data.size
        inds_data = inds_data[m]
        inds_repeats = np.where(m)[0]
        m = offset_num_data['repeat_index'][inds_data] == repeat_indexes[inds_repeats]
        inds_data = inds_data[m]
        return offset_num_data[inds_data]
        

    def repeat_indexes_from_offset(self, offset) :
        offset_num_data = self.num_data_from_offset(offset)
        return offset_num_data['repeat_index']

    def num_data_from_offset_range(self, offset_start, offset_bound) :
        m = self.num_data['alu_offset'] >= offset_start
        m = np.logical_and(m, self.num_data['alu_offset'] < offset_bound)
        return self.num_data[m]

    def repeats_from_indexes(self, repeat_indexes) :
        inds_rd = self.repeats['index'].searchsorted(repeat_indexes)
        return self.repeats[inds_rd]
    
    def num_count(self, num, offset=None) :
        m = self.num_data['num_16'] == num
        if offset is not None :
            mo = self.num_data['alu_offset'] == offset
            m = np.logical_and(m, mo)
        return m.sum()
    
    def num_data_from_num(self, num_16, offset=None) :
        m = self.num_data['num_16'] == num_16
        if offset is not None :
            mo = self.num_data['alu_offset'] == offset
            m = np.logical_and(m, mo)
        return self.num_data[m]

    def offset_counts_from_num(self, num_16) :
        num_data = self.num_data_from_num(num_16)
        offsets, counts = np.unique(num_data['alu_offset'], return_counts=True)
        data = np.zeros(offsets.size, dtype=self.offset_count_dtype)
        data['alu_offset'] = offsets
        data['count'] = counts
        return data
    
    def num_counts_from_offset(self, offset) :
        num_data = self.num_data_from_offset(offset)
        nums, counts = np.unique(num_data['num_16'], return_counts=True)
        data = np.zeros(nums.size, self.num_count_dtype)
        data['num_16'] = nums
        data['count'] = counts
        return data
        
    def num_counts_from_offset_range(self, offset_start, offset_bound) :
        num_data = self.num_data_from_offset_range(offset_start, offset_bound)
        nums, counts = np.unique(num_data['num_16'], return_counts=True)
        data = np.zeros(nums.size, self.num_count_dtype)
        data['num_16'] = nums
        data['count'] = counts
        return data

        

    def cg_repeat_indexes_from_offset(self, offset) :
        num_data = self.num_data_from_offset(offset)
        m_cg = n16d.cg_high_mask(num_data['num_16'])
        return num_data['repeat_index'][m_cg]

    def or_cg_tg_ca_repeat_indexes_from_offset(self, offset) :
        num_data = self.num_data_from_offset(offset)
        nums = num_data['num_16']
        m = n16d.cg_high_mask(nums)
        m = np.logical_or(m, n16d.tg_high_mask(nums))
        m = np.logical_or(m, n16d.ca_high_mask(nums))
        return num_data['repeat_index'][m]

    def obj_from_repeat_indexes(self, repeat_indexes) :
        num_data = self.num_data_from_repeat_indexes(repeat_indexes)
        inds_repeats = self.repeats['index'].searchsorted(repeat_indexes)
        repeats = self.repeats[inds_repeats]
        return repeat_num_data_cls(num_data, repeats)
        
    def seq_str_from_repeat_index(self, repeat_index) :
        bound = repeat_index + 1
        inds = self.num_data['repeat_index'].searchsorted([repeat_index, bound])
        num_data = self.num_data[inds[0]:inds[1]]
        seq_str = n16d.dna_str_from_num_16_array(num_data['num_16'])
        return seq_str

    def disply_num_counts_from_offset(self, offset, min_count=100) :
        num_counts = self.num_counts_from_offset(offset)
        m = num_counts['count'] >= min_count
        data_size = m.sum()
        oda = np.zeros(data_size, dtype=self.num_count_letters_dtype)
        oda['num_16'] = num_counts['num_16'][m]
        oda['count'] = num_counts['count'][m]
        oda['letters'] = n16d.num_16_strs(oda['num_16'])
        oda.sort(order='count')
        oda = oda[::-1]
        return oda
        
class strand_alu_data_base_cls(object) :
    local_folder = 'grch37_hg19_alu_data'
    data_folder = os.path.join(mod_dir, local_folder)
    
    def __init__(self) :
        self.read_data()
        self.num_counts = None
        self.offset_stats = None
        self.cg_counts_obj = None

    def read_data(self) :
        file_path = os.path.join(self.data_folder, self.file_name)
        h5 = tb.open_file(file_path, 'r')
        num_data_table = getattr(h5.root, self.num_data_table_name)
        self.num_data = num_data_table[:]
        repeat_table = getattr(h5.root, self.repeat_table_name)
        self.repeats = repeat_table[:]
        self.repeat_count = self.repeats.size
        h5.close()
        self.num_data_obj = repeat_num_data_cls(self.num_data, self.repeats)

    def read_counts_and_stats(self) :
        if self.num_counts is None :
            self.read_num_counts()
        if self.offset_stats is None :
            self.read_offset_stats()
            
    def data_all_repeat_indexes(self) :
        return self.repeats['index']

    def num_data_from_repeat_indexes(self, repeat_indexes) :
        bounds = repeat_indexes + 1
        inds_starts = self.num_data['repeat_index'].searchsorted(repeat_indexes)
        inds_bounds = self.num_data['repeat_index'].searchsorted(bounds)
        out_data = []
        for ind_start, ind_bound in zip(inds_starts, inds_bounds) :
            out_data.append(self.num_data[ind_start:ind_bound])
        out_data = np.concatenate(out_data)
        return out_data            
            
    def repeat_indexes_from_num(self, num, offset=None) :
        m = self.num_data['num_16'] == num
        if offset is not None :
            mo = self.num_data['alu_offset'] == offset
            m = np.logical_and(m, mo)
        repeat_indexes = self.num_data['repeat_index'][m]
        if offset is None :
            repeat_indexes = np.unique(repeat_indexes)
        return repeat_indexes

    def repeat_indexes_from_offset_nums(self, offset_nums) :
        ris = None
        for offset, num in offset_nums :
            oris = self.repeat_indexes_from_num(num, offset=offset)
            if ris is None :
                ris = oris
            else :
                ris = np.intersect1d(ris, oris)
        return ris

    def num_data_from_offset(self, offset) :
        m = self.num_data['alu_offset'] == offset
        return self.num_data[m]
    
    def repeats_from_indexes(self, repeat_indexes) :
        inds_rd = self.repeats['index'].searchsorted(repeat_indexes)
        return self.repeats[inds_rd]
    
    def num_count(self, num, offset=None) :
        m = self.num_data['num_16'] == num
        if offset is not None :
            mo = self.num_data['alu_offset'] == offset
            m = np.logical_and(m, mo)
        return m.sum()
    
    def num_data_from_num(self, num_16, offset=None) :
        m = self.num_data['num_16'] == num_16
        if offset is not None :
            mo = self.num_data['alu_offset'] == offset
            m = np.logical_and(m, mo)
        return self.num_data[m]
        
    def obj_from_repeat_indexes(self, repeat_indexes) :
        num_data = self.num_data_from_repeat_indexes(repeat_indexes)
        inds_repeats = self.repeats['index'].searchsorted(repeat_indexes)
        repeats = self.repeats[inds_repeats]
        return repeat_num_data_cls(num_data, repeats)
        
class pos_alu_data_cls(strand_alu_data_base_cls) :
    file_name = 'pos_alu_dna_num_16.h5'
    num_data_table_name = 'pos_alu_dna_num_16_data'
    repeat_table_name = 'alu_plus_strand_repeats'

    def __init__(self) :
        strand_alu_data_base_cls.__init__(self)
        
    def read_num_counts(self) :
        num_counts_obj = pos_alu_num_counts()
        self.num_counts = num_counts_obj.num_counts

    def read_offset_stats(self) :
        offset_stats_obj = pos_alu_offset_stats_cls()
        self.offset_stats = offset_stats_obj.offset_stats
        
    def read_cg_counts_obj(self) :
        if self.cg_counts_obj is None :
            self.cg_counts_obj = pos_alu_cg_cls()
        
class neg_alu_data_cls(strand_alu_data_base_cls) :
    file_name = 'neg_alu_dna_num_16.h5'
    num_data_table_name = 'neg_alu_dna_num_16_data'
    repeat_table_name = 'alu_minus_strand_repeats'

    def __init__(self) :
        strand_alu_data_base_cls.__init__(self)
        
    def read_num_counts(self) :
        num_counts_obj = neg_alu_num_counts()
        self.num_counts = num_counts_obj.num_counts
        
    def read_offset_stats(self) :
        offset_stats_obj = neg_alu_offset_stats_cls()
        self.offset_stats = offset_stats_obj.offset_stats
        
    def read_cg_counts_obj(self) :
        if self.cg_counts_obj is None :
            self.cg_counts_obj = neg_alu_cg_cls()
        
        
class alu_num_counts_cls(object) :
    local_folder = 'grch37_hg19_alu_data'
    data_folder = os.path.join(mod_dir, local_folder)
    count_file_name = 'merged_num_16_counts.h5'
    file_path = os.path.join(data_folder, count_file_name)
    count_table_name = 'merged_num_16_counts'
    
    def __init__(self) :
        self.read_data()
    
    def read_data(self) :
        h5 = tb.open_file(self.file_path, 'r')
        count_table = getattr(h5.root, self.count_table_name)
        self.num_counts = count_table[:]
        h5.close()
        
    def counts_from_nums(self, nums) :
        inds_counts = self.num_counts['num_16'].searchsorted(nums)
        return self.counts[inds_counts]
    
class pos_alu_num_counts(alu_num_counts_cls) :
    count_table_name = 'pos_num_16_counts'
    
    def __init__(self) :
        alu_num_counts_cls.__init__(self)


class neg_alu_num_counts(alu_num_counts_cls) :
    count_table_name = 'neg_num_16_counts'
    
    def __init__(self) :
        alu_num_counts_cls.__init__(self)


class alu_strand_offset_stats_base_cls(object) :
    local_folder = 'grch37_hg19_alu_data'
    data_folder = os.path.join(mod_dir, local_folder)
    file_name = 'alu_strand_offset_data_stats.h5'
    file_path = os.path.join(data_folder, file_name)
    
    def __init__(self) :
        self.read_data()
        
    def read_data(self) :
        h5 = tb.open_file(self.file_path, 'r')
        stats_table = getattr(h5.root, self.stats_table_name)
        self.offset_stats = stats_table[:]
        h5.close()
        
class pos_alu_offset_stats_cls(alu_strand_offset_stats_base_cls) :
    stats_table_name = 'pos_alu_offset_data_stats'
    
    def __init__(self) :
        alu_strand_offset_stats_base_cls.__init__(self)
        
class neg_alu_offset_stats_cls(alu_strand_offset_stats_base_cls) :
    stats_table_name = 'neg_alu_offset_data_stats'
    
    def __init__(self) :
        alu_strand_offset_stats_base_cls.__init__(self)
        
        
class alu_strand_cg_base_cls(object) :        
    local_folder = 'grch37_hg19_alu_data'
    data_folder = os.path.join(mod_dir, local_folder)
    cg_count_dtype = np.dtype([('cg_count', np.uint32), ('count', np.uint32)])
    
    def __init__(self) :
        self.read_data()
        
    def read_data(self) :
        file_path = os.path.join(self.data_folder, self.file_name)
        h5 = tb.open_file(file_path, 'r')
        repeat_table = getattr(h5.root, self.alu_repeat_table_name)
        self.repeats = repeat_table[:]
        repeat_index_cg_offset_table = getattr(h5.root, self.repeat_index_cg_offset_table_name)
        self.repeat_indexes_cg_offsets = repeat_index_cg_offset_table[:]
        offset_cg_count_table = getattr(h5.root, self.offset_cg_count_table_name)
        self.offset_cg_counts = offset_cg_count_table[:]
        dist_repeat_cg_counts_table = getattr(h5.root, self.dist_repeat_cg_counts_table_name)
        self.dist_repeat_cg_counts = dist_repeat_cg_counts_table[:]
        h5.close()
        
    def repeats_from_indexes(self, repeat_indexes) :
        inds_rd = self.repeats['index'].searchsorted(repeat_indexes)
        return self.repeats[inds_rd]

    def repeat_cg_counts(self, name=None) :
        if name is None :
            repeats = self.repeats
        else :
            m = self.repeats['repeat_name'] == name
            repeats = self.repeats[m]
        cg_counts, counts = np.unique(repeats['cg_count'], return_counts=True)
        od = np.zeros(cg_counts.size, self.cg_count_dtype)
        od['cg_count'] = cg_counts
        od['count'] = counts
        return od
    
    def repeats_by_cg_count(self, min_cg_count=None, max_cg_count=None) :
        if min_cg_count :
            repeat_mask = self.repeats['cg_count'] >= min_cg_count
        else :
            repeat_mask = None
        if max_cg_count :
            max_repeat_mask = self.repeats['cg_count'] < max_cg_count
            if repeat_mask :
                repeat_mask = np.logical_and(repeat_mask, max_repeat_mask)
            else :
                repeat_mask = max_repeat_mask
        return self.repeats['index'][repeat_mask]

    def repeat_indexes_from_offset(self, offset) :
        m = self.repeat_indexes_cg_offsets['alu_offset'] == offset
        return self.repeat_indexes_cg_offsets['repeat_index'][m]
                
class pos_alu_cg_cls(alu_strand_cg_base_cls) :
    file_name = 'pos_strand_repeat_cg_data.h5'    
    alu_repeat_table_name = 'alu_plus_strand_repeats'
    repeat_index_cg_offset_table_name = 'pos_repeat_cg_offsets'
    offset_cg_count_table_name = 'pos_offset_cg_counts'
    dist_repeat_cg_counts_table_name = 'pos_repeat_cg_count_distribution'
    
    def __init__(self) :
        alu_strand_cg_base_cls.__init__(self)
        
        
class neg_alu_cg_cls(alu_strand_cg_base_cls) :
    file_name = 'neg_strand_repeat_cg_data.h5'    
    alu_repeat_table_name = 'alu_minus_strand_repeats'
    repeat_index_cg_offset_table_name = 'neg_repeat_cg_offsets'
    offset_cg_count_table_name = 'neg_offset_cg_counts'
    dist_repeat_cg_counts_table_name = 'neg_repeat_cg_count_distribution'
    
    def __init__(self) :
        alu_strand_cg_base_cls.__init__(self)
        
        
class alu_clusters_base_cls(object) :
    local_folder = 'grch37_hg19_alu_data'
    data_folder = os.path.join(mod_dir, local_folder)

    def __init__(self) :
        self.read_data()
    
    def read_data(self) :
        file_path = os.path.join(self.data_folder, self.cluster_file_name)
        h5 = tb.open_file(file_path, 'r')
        cluster_table = getattr(h5.root, self.cluster_table_name)
        self.cluster_num_counts = cluster_table[:]
        poly_table = getattr(h5.root, self.poly_table_name)
        self.poly_num_counts = poly_table[:]
        h5.close()

        
class pos_alu_clusters_cls(alu_clusters_base_cls) :
    cluster_file_name = 'pos_alu_num_16_clusters.h5'
    cluster_table_name = 'pos_alu_num_16_clusters'
    poly_table_name = 'pos_alu_repeat_poly_num_16s'
    
    def __init__(self) :
        alu_clusters_base_cls.__init__(self)
        
        
class neg_alu_clusters_cls(alu_clusters_base_cls) :
    cluster_file_name = 'neg_alu_num_16_clusters.h5'
    cluster_table_name = 'neg_alu_num_16_clusters'
    poly_table_name = 'neg_alu_repeat_poly_num_16s'
    
    def __init__(self) :
        alu_clusters_base_cls.__init__(self)
        


class alu_offset_num_counts_base_cls(object) :        
    local_folder = 'grch37_hg19_alu_data'
    data_folder = os.path.join(mod_dir, local_folder)
        
    def __init__(self) :
        self.read_data()

    def read_data(self) :        
        file_path = os.path.join(self.data_folder, self.offset_num_counts_file_name)
        h5 = tb.open_file(file_path, 'r')
        offset_num_counts_table = getattr(h5.root, self.offset_num_counts_table_name)
        self.offset_num_counts = offset_num_counts_table[:]
        summary_num_counts_table = getattr(h5.root, self.summary_num_counts_table_name)
        self.num_offset_summary_data = summary_num_counts_table[:]

    def offsets_from_num(self, num_16) :
        bound = num_16 + 1
        inds = self.offset_num_counts['num_16'].searchsorted([num_16, bound])
        nos = self.offset_num_counts[inds[0]:inds[1]].copy()
        nos.sort(order='alu_offset')
        return nos

        
class pos_offset_num_counts_cls(alu_offset_num_counts_base_cls) :
    offset_num_counts_file_name = 'pos_alu_offset_num_counts.h5'
    offset_num_counts_table_name = 'pos_alu_offset_num_counts'    
    summary_num_counts_table_name = 'pos_alu_num_offset_count_summaries'    
    
    def __init__(self) :
        alu_offset_num_counts_base_cls.__init__(self)
        
        
class neg_offset_num_counts_cls(alu_offset_num_counts_base_cls) :
    offset_num_counts_file_name = 'neg_alu_offset_num_counts.h5'
    offset_num_counts_table_name = 'neg_alu_offset_num_counts'    
    summary_num_counts_table_name = 'neg_alu_num_offset_count_summaries'    
    
    def __init__(self) :
        alu_offset_num_counts_base_cls.__init__(self)
        
        


        
    

