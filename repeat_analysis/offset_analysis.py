# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d
import genomes_dnj_2.repeat_analysis.parital_sequences_analysis as psa



class offset_anal_cls(object) :
    num_count_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32)])
    offset_dtype = np.dtype([('offset', np.uint32), ('num_16', np.uint32), ('letters', 'S16'),
                             ('offset_max_num_count', np.uint32), ('offset_total_num_count', np.uint32),
                             ('offset_unique_nums', np.uint32), ('anal_max_num_count', np.uint32),
                             ('strand_max_num_offset_count', np.uint32), ('strand_max_num_total_count', np.uint32),
                             ('strand_max_num_repeat_count', np.uint32), ('rel_expr_max_num_offset_count', np.float32),
                             ('rel_expr_max_num_total_count', np.float32)])
    
    def __init__(self, num_data, data_obj) :
        self.num_data = num_data
        self.data_obj = data_obj
        #self.strand = strand
        self.data_obj.read_counts_and_stats()
        
    def count_num_16s(self) :
        nums, counts = np.unique(self.num_data['num_16'], return_counts=True)
        self.anal_num_counts = np.zeros(nums.size, dtype=self.num_count_dtype)
        self.anal_num_counts['num_16'] = nums
        self.anal_num_counts['count'] = counts

    def process_offset(self, offset, offset_num_data) :
        offset_nums, offset_counts = np.unique(offset_num_data['num_16'], return_counts=True)
        ami = offset_counts.argmax()
        num_16 = offset_nums[ami]
        letters = ''
        offset_max_num_count = offset_counts[ami]
        offset_total_num_count = offset_counts.sum()
        offset_unique_nums = offset_nums.size
        anal_max_num_count = 0
        strand_max_num_offset_count = self.data_obj.num_count(num_16, offset=offset)
        strand_max_num_total_count = 0
        strand_max_num_repeat_count = 0
        rel_expr_max_num_offset_count = 0.0
        rel_expr_max_num_total_count = 0.0
        item = (offset, num_16, letters, offset_max_num_count, offset_total_num_count, offset_unique_nums,
                anal_max_num_count, strand_max_num_offset_count, strand_max_num_total_count, strand_max_num_repeat_count,
                rel_expr_max_num_offset_count, rel_expr_max_num_total_count)
        return item
    
    def process_offsets(self) :
        self.num_data.sort(order=['alu_offset', 'num_16'])
        offsets, starts, counts = np.unique(self.num_data['alu_offset'], return_index=True, return_counts=True)
        out_data = []
        for offset, start, count in zip(offsets, starts, counts) :
            print('processing', offset)
            bound = start + count
            offset_num_data = self.num_data[start:bound]
            data_item = self.process_offset(offset, offset_num_data)
            out_data.append(data_item)
        self.offset_num_counts = np.array(out_data, dtype=self.offset_dtype)

    def calculate_rel_expr(self) :
        ac = self.anal_num_counts
        ac_nums = ac['num_16']
        onc = self.offset_num_counts
        onc_nums = onc['num_16']
        aci = ac_nums.searchsorted(onc_nums)
        onc['anal_max_num_count'] = ac['count'][aci]
        dnc = self.data_obj.num_counts
        dnc_nums = dnc['num_16']
        dci = dnc_nums.searchsorted(onc_nums)
        onc['strand_max_num_total_count'] = dnc['total_count'][dci]
        onc['strand_max_num_repeat_count'] = dnc['repeat_count'][dci]
        data_repeat_count = float(self.data_obj.repeat_count)
        anal_offset_repeat_counts = onc['offset_total_num_count'].astype(np.float32)
        anal_offset_max_num_counts = onc['offset_max_num_count'].astype(np.float32)
        anal_offset_max_num_total_count = onc['anal_max_num_count'].astype(np.float32)
        data_offset_max_num_repeat_counts = onc['strand_max_num_repeat_count'].astype(np.float32)
        data_offset_max_num_offset_counts = onc['strand_max_num_offset_count'].astype(np.float32)
        offset_pred_factors = anal_offset_repeat_counts/data_repeat_count
        pred_offset_max_num_offset_counts = offset_pred_factors*data_offset_max_num_offset_counts
        pred_offset_max_num_total_counts = offset_pred_factors*data_offset_max_num_repeat_counts
        onc['rel_expr_max_num_offset_count'] = anal_offset_max_num_counts/pred_offset_max_num_offset_counts
        onc['rel_expr_max_num_total_count'] = anal_offset_max_num_total_count/pred_offset_max_num_total_counts
        onc['letters'] = n16d.num_16_strs(onc_nums)

    def do_work(self) :
        self.count_num_16s()
        self.process_offsets()
        self.calculate_rel_expr()
        
'''
Probably want to make the class initialization just depend on repeat indexes
because I want a number of ways to choose the focus set of repeats
'''        
        
class offset_num_repeats_anal_cls(object) :
    data_dtype = np.dtype([('offset', np.uint16), ('offset_max_num', np.uint32), ('letters', 'S16'), 
                           ('offset_max_num_count', np.uint32), ('offset_total_count', np.uint32), 
                           ('offset_total_nums', np.uint32), ('all_offset_max_num_count', np.uint32),
                           ('stats_offset_repeat_count', np.uint32), ('rel_expr_offset_max_num', np.float32), ('stats_max_num', np.uint32), 
                           ('stats_max_num_count', np.uint32), ('offset_stats_max_num_count', np.uint32),
                           ('rel_expr_stats_max_num', np.float32), ('offset_cg_count', np.uint32),
                           ('stats_offset_cg_count', np.uint32), ('rel_expr_cg_count', np.float32),
                           ('offset_cg_tg_ca_count', np.uint32), ('stats_offset_cg_tg_ca_count', np.uint32), ('rel_expr_cg_tg_ca_count', np.float32)])
    rel_num_count_dtype = np.dtype([('num_16', np.uint32), ('repeat_offset_num_count', np.uint32), 
                                    ('data_offset_num_count', np.uint32), ('rel_expr', np.float32), ('letters', 'S16')])
    rel_offset_count_dtype = np.dtype([('alu_offset', np.uint32), ('repeat_num_offset_count', np.uint32), 
                                    ('data_num_offset_count', np.uint32), ('rel_expr', np.float32)])
    offset_cg_dtype = np.dtype([('offset', np.uint16), ('offset_cg_count', np.uint32), 
                                ('stats_offset_cg_count', np.uint32), ('rel_expr_cg_count', np.float32)])
    offset_cg_derived_dtype = np.dtype([('offset', np.uint16), ('offset_cg_tg_ca_count', np.uint32), 
                                ('stats_offset_cg_tg_ca_count', np.uint32), ('rel_expr_cg_tg_ca_count', np.float32)])

    @classmethod
    def from_num_offset(cls, num_16, offset, data_obj) :
        repeat_indexes = data_obj.repeat_indexes_from_num(num_16, offset)
        return cls(repeat_indexes, data_obj)
    
    @classmethod
    def from_offset_num_pairs(cls, offset_nums, data_obj) :
        repeat_indexes = None
        for offset, num in offset_nums :
            n_o_repeat_indexes = data_obj.repeat_indexes_from_num(num, offset=offset)
            if repeat_indexes is None :
                repeat_indexes = n_o_repeat_indexes
            else :
                repeat_indexes = np.intersect1d(repeat_indexes, n_o_repeat_indexes)
        return cls(repeat_indexes, data_obj)
    
    @classmethod
    def from_yes_no_offset_num_pairs(cls, yes_offset_nums, no_offset_nums, data_obj) :
        repeat_indexes = None
        for offset, num in yes_offset_nums :
            n_o_repeat_indexes = data_obj.repeat_indexes_from_num(num, offset=offset)
            if repeat_indexes is None :
                repeat_indexes = n_o_repeat_indexes
            else :
                repeat_indexes = np.intersect1d(repeat_indexes, n_o_repeat_indexes)
                
        for offset, num in no_offset_nums :
            n_o_repeat_indexes = data_obj.repeat_indexes_from_num(num, offset=offset)
            repeat_indexes = np.setdiff1d(repeat_indexes, n_o_repeat_indexes)
        
        return cls(repeat_indexes, data_obj)

    @classmethod
    def from_num(cls, num_16, data_obj) :
        repeat_indexes = data_obj.repeat_indexes_from_num(num_16)
        return cls(repeat_indexes, data_obj)
    
    @classmethod
    def from_num_and_repeat_indexes(cls, num_16, repeat_indexes, data_obj) :
        num_repeat_indexes = data_obj.repeat_indexes_from_num(num_16)
        repeat_indexes = np.intersect1d(repeat_indexes, num_repeat_indexes)
        return cls(repeat_indexes, data_obj)

    @classmethod
    def from_alu_name(cls, alu_name, data_obj) :
        m = data_obj.repeats['repeat_name'] == alu_name
        repeat_indexes = data_obj.repeats['index'][m]
        return cls(repeat_indexes, data_obj)
    
    @classmethod
    def from_alu_names(cls, alu_names, data_obj) :
        repeat_indexes = None
        for name in alu_names :
            m = data_obj.repeats['repeat_name'] == name
            name_repeat_indexes = data_obj.repeats['index'][m]
            if repeat_indexes is None :
                repeat_indexes = name_repeat_indexes
            else :
                repeat_indexes = np.union1d(repeat_indexes, name_repeat_indexes)
        return cls(repeat_indexes, data_obj)

    @classmethod
    def from_cg_offset(cls, offset, data_obj, in_repeat_indexes=None) :
        num_data = data_obj.num_data_obj.num_data_from_offset(offset)
        m = n16d.cg_high_mask(num_data['num_16'])
        repeat_indexes = num_data['repeat_index'][m]
        if in_repeat_indexes is not None :
            repeat_indexes = np.intersect1d(repeat_indexes, in_repeat_indexes)
        return cls(repeat_indexes, data_obj) 
            
    @classmethod
    def from_and_cg_offsets(cls, data_obj, yes_offsets=None, no_offsets=None, in_repeat_indexes=None, ) :
        repeat_indexes = None
        if yes_offsets is None :
            yes_offsets = []
            if in_repeat_indexes is None :
                repeat_indexes = data_obj.repeats['index']
            else :
                repeat_indexes = in_repeat_indexes
        for offset in yes_offsets :
            num_data = data_obj.num_data_obj.num_data_from_offset(offset)
            m = n16d.cg_high_mask(num_data['num_16'])        
            offset_repeat_indexes = num_data['repeat_index'][m]
            if repeat_indexes is None :
                repeat_indexes = offset_repeat_indexes
            else :
                repeat_indexes = np.intersect1d(repeat_indexes, offset_repeat_indexes)
        if no_offsets is None :
            no_offsets = []
        for offset in no_offsets :
            num_data = data_obj.num_data_obj.num_data_from_offset(offset)
            m = n16d.cg_high_mask(num_data['num_16'])        
            offset_repeat_indexes = num_data['repeat_index'][m]
            repeat_indexes = np.setdiff1d(repeat_indexes, offset_repeat_indexes)                            
        if in_repeat_indexes is not None :
            repeat_indexes = np.intersect1d(repeat_indexes, in_repeat_indexes)
        return cls(repeat_indexes, data_obj) 

    @classmethod
    def from_and_cg_tg_ca_offsets(cls, data_obj, yes_offsets=None, no_offsets=None, 
                                  in_repeat_indexes=None, alu_name=None, add_gg=False ) :
        repeat_indexes = None
        if yes_offsets is None :
            yes_offsets = []
            if in_repeat_indexes is None :
                repeat_indexes = data_obj.repeats['index']
            else :
                repeat_indexes = in_repeat_indexes
        for offset in yes_offsets :
            num_data = data_obj.num_data_obj.num_data_from_offset(offset)
            if add_gg :
                m = n16d.or_cg_tg_ca_gg_high_mask(num_data['num_16'])
            else :
                m = n16d.or_cg_tg_ca_high_mask(num_data['num_16'])
            offset_repeat_indexes = num_data['repeat_index'][m]
            if repeat_indexes is None :
                repeat_indexes = offset_repeat_indexes
            else :
                repeat_indexes = np.intersect1d(repeat_indexes, offset_repeat_indexes)
        if no_offsets is None :
            no_offsets = []
        for offset in no_offsets :
            num_data = data_obj.num_data_obj.num_data_from_offset(offset)
            if add_gg :
                m = n16d.or_cg_tg_ca_gg_high_mask(num_data['num_16'])
            else :   
                m = n16d.or_cg_tg_ca_high_mask(num_data['num_16'])        
            offset_repeat_indexes = num_data['repeat_index'][m]
            repeat_indexes = np.setdiff1d(repeat_indexes, offset_repeat_indexes)                            
        if in_repeat_indexes is not None :
            repeat_indexes = np.intersect1d(repeat_indexes, in_repeat_indexes)
        if alu_name is not None :
            m = data_obj.repeats['repeat_name'] == alu_name
            name_repeat_indexes = data_obj.repeats['index'][m]
            repeat_indexes = np.intersect1d(repeat_indexes, name_repeat_indexes)
            
        return cls(repeat_indexes, data_obj) 



    def __init__(self, repeat_indexes, data_obj) :
        self.repeat_indexes = repeat_indexes
        self.data_obj = data_obj
        self.num_data_obj = self.data_obj.obj_from_repeat_indexes(self.repeat_indexes)
        self.data_obj.read_counts_and_stats()
        self.offset_stats = self.data_obj.offset_stats
        self.data_obj.read_cg_counts_obj()
        self.cg_stats = self.data_obj.cg_counts_obj.offset_cg_counts
        self.pred_factor = float(self.repeat_indexes.size)/float(self.data_obj.repeat_count)
        
    def process_offset(self, offset) :
         offset_num_data = self.num_data_obj.num_data_from_offset(offset)
         offset_total_count = offset_num_data.size
         nums, counts = np.unique(offset_num_data['num_16'], return_counts=True)
         offset_total_nums = nums.size
         if nums.size == 0 :
             return None
         imax = counts.argmax()
         offset_max_num = nums[imax]
         offset_max_num_count = counts[imax]
         offset_stats_index = self.offset_stats['alu_offset'].searchsorted(offset)
         offset_stats = self.offset_stats[offset_stats_index]
         offset_cg_stats_index = self.cg_stats['alu_offset'].searchsorted(offset)
         offset_cg_stats = self.cg_stats[offset_cg_stats_index]
         stats_offset_repeat_count = offset_stats['total_offset_count']
         data_offset_max_num_data = self.data_obj.num_data_from_num(offset_max_num, offset=offset)         
         all_offset_max_num_count = data_offset_max_num_data.size
         offset_pred_factor = float(offset_total_count)/float(stats_offset_repeat_count)
         pred_offset_max_num_count = offset_pred_factor*all_offset_max_num_count
         rel_expr_offset_max_num = float(offset_max_num_count)/pred_offset_max_num_count                  
         stats_max_num = offset_stats['max_num_16']
         ind_stats_max_num = nums.searchsorted(stats_max_num)
         offset_stats_max_num_count = 0
         if ((ind_stats_max_num < nums.size) and
             (stats_max_num == nums[ind_stats_max_num])) :
                 offset_stats_max_num_count = counts[ind_stats_max_num]             
         stats_max_num_count = offset_stats['max_num_offset_count']                  
         stats_pred_factor = float(stats_max_num_count)/float(stats_offset_repeat_count)
         pred_offset_stats_max_num_count = stats_pred_factor*float(offset_total_count)
         rel_expr_stats_max_num = float(offset_stats_max_num_count)/pred_offset_stats_max_num_count
         offset_cg_mask = n16d.cg_high_mask(nums)
         offset_cg_counts = counts[offset_cg_mask]
         offset_cg_count = offset_cg_counts.sum()
         offset_tg_mask = n16d.tg_high_mask(nums)
         offset_tg_counts = counts[offset_tg_mask]
         offset_tg_count = offset_tg_counts.sum()
         offset_ca_mask = n16d.ca_high_mask(nums)
         offset_ca_counts = counts[offset_ca_mask]
         offset_ca_count = offset_ca_counts.sum()
         stats_offset_cg_count = offset_cg_stats['cg_count']
         stats_offset_tg_count = offset_cg_stats['tg_count']
         stats_offset_ca_count = offset_cg_stats['ca_count']
         if stats_offset_cg_count > 0 :
             pred_cg_count = offset_pred_factor*float(stats_offset_cg_count)
             rel_expr_cg_count = float(offset_cg_count)/pred_cg_count
         else :
             rel_expr_cg_count = 0.0
         offset_all_count = offset_cg_count + offset_tg_count + offset_ca_count
         stats_offset_all_count = stats_offset_cg_count + stats_offset_tg_count + stats_offset_ca_count
         if stats_offset_all_count > 0 :
             pred_all_count = offset_pred_factor*float(stats_offset_all_count)
             rel_expr_all_count = float(offset_all_count)/pred_all_count
         data = (offset, offset_max_num, '', offset_max_num_count, offset_total_count, offset_total_nums,
                 all_offset_max_num_count, stats_offset_repeat_count, rel_expr_offset_max_num, stats_max_num,
                 stats_max_num_count, offset_stats_max_num_count, rel_expr_stats_max_num,
                 offset_cg_count, stats_offset_cg_count, rel_expr_cg_count,
                 offset_all_count, stats_offset_all_count, rel_expr_all_count)
         return data
     
    def process_offsets(self) :
        out_data = []
        for offset in range(1, 300) :
            offset_data = self.process_offset(offset)
            if offset_data is not None :
                out_data.append(offset_data)
        self.offset_counts = np.array(out_data, dtype=self.data_dtype)
        self.offset_counts['letters'] = n16d.num_16_strs(self.offset_counts['offset_max_num'])

    def most_likely_sequence_string(self) :
        self.sequence_string = n16d.dna_str_from_num_16_array(self.offset_counts['offset_max_num'])
    
    def ttacg_repeat_indexes(self, offset) :
        offset_num_data = self.num_data_obj.num_data_from_offset(offset)
        mdm_obj = psa.masked_num_match_cls(offset_num_data['num_16'], 5)
        m = mdm_obj.match_num(psa.masked_ttagc)
        matched_ris = offset_num_data['repeat_index'][m]
        return matched_ris

    def relative_num_counts(self, offset, min_data_count=100) :
        repeat_num_counts = self.num_data_obj.num_counts_from_offset(offset)
        data_num_counts = self.data_obj.num_data_obj.num_counts_from_offset(offset)
        m = data_num_counts['count'] >= min_data_count
        count_data_size = m.sum()
        oda = np.zeros(count_data_size, dtype=self.rel_num_count_dtype)
        oda['num_16'] = data_num_counts['num_16'][m]
        oda['data_offset_num_count'] = data_num_counts['count'][m]
        inds = oda['num_16'].searchsorted(repeat_num_counts['num_16'])
        m = inds < oda.size
        m[m] = oda['num_16'][inds[m]] == repeat_num_counts['num_16'][m]
        oda['repeat_offset_num_count'][inds[m]] = repeat_num_counts['count'][m]
        pred_repeat_offset_num_counts = self.pred_factor*oda['data_offset_num_count'].astype(np.float32)
        oda['rel_expr'] = oda['repeat_offset_num_count'].astype(np.float32)/pred_repeat_offset_num_counts
        oda['letters'] = n16d.num_16_strs(oda['num_16'])
        oda.sort(order='data_offset_num_count')
        oda = oda[::-1]
        return oda

    def repeat_sorted_rel_num_counts(self, offset, min_data_count=100) :
        rnc = self.relative_num_counts(offset, min_data_count)
        rnc.sort(order=['repeat_offset_num_count', 'data_offset_num_count'])
        rnc = rnc[::-1]
        return rnc

    def relative_offset_counts(self, num_16) :
        repeat_offset_counts = self.num_data_obj.offset_counts_from_num(num_16)
        data_offset_counts = self.data_obj.num_data_obj.offset_counts_from_num(num_16)
        oda = np.zeros(data_offset_counts.size, dtype=self.rel_offset_count_dtype)
        oda['alu_offset'] = data_offset_counts['alu_offset']
        oda['data_num_offset_count'] = data_offset_counts['count']
        inds = oda['alu_offset'].searchsorted(repeat_offset_counts['alu_offset'])
        oda['repeat_num_offset_count'][inds] = repeat_offset_counts['count']
        pred_repeat_num_offset_counts = self.pred_factor*oda['data_num_offset_count'].astype(np.float32)
        oda['rel_expr'] = oda['repeat_num_offset_count'].astype(np.float32)/pred_repeat_num_offset_counts
        return oda

    
    def output_masked_counts(self, m) :
        oc = self.offset_counts
        od = np.zeros(m.sum(), dtype=self.offset_cg_dtype)
        od['offset'] = oc['offset'][m]
        od['offset_cg_count'] = oc['offset_cg_count'][m]
        od['stats_offset_cg_count'] = oc['stats_offset_cg_count'][m]
        od['rel_expr_cg_count'] = oc['rel_expr_cg_count'][m]
        return od
        
    def output_derived_masked_counts(self, m) :
        oc = self.offset_counts
        od = np.zeros(m.sum(), dtype=self.offset_cg_derived_dtype)
        od['offset'] = oc['offset'][m]
        od['offset_cg_tg_ca_count'] = oc['offset_cg_tg_ca_count'][m]
        od['stats_offset_cg_tg_ca_count'] = oc['stats_offset_cg_tg_ca_count'][m]
        od['rel_expr_cg_tg_ca_count'] = oc['rel_expr_cg_tg_ca_count'][m]
        return od

    def high_cg_derived_counts(self, min_count) :
        oc = self.offset_counts        
        m = oc['offset_cg_tg_ca_count'] >= min_count
        #return self.output_derived_masked_counts(m)
        return self.offset_counts[m]
    
    def high_cg_derived_rel_expr(self, min_rel_expr) :
        oc = self.offset_counts        
        m = oc['rel_expr_cg_tg_ca_count'] >= min_rel_expr
        #return self.output_derived_masked_counts(m)
        return self.offset_counts[m]

    def high_cg_counts(self, min_count) :
        oc = self.offset_counts        
        m = oc['offset_cg_count'] >= min_count
        #return self.output_masked_counts(m)
        return self.offset_counts[m]
    
    def high_cg_rel_expr(self, min_rel_expr) :
        oc = self.offset_counts        
        m = oc['rel_expr_cg_count'] >= min_rel_expr
        #return self.output_masked_counts(m)
        return self.offset_counts[m]
        


        
class num_offsets_anal_cls(object) :
    offset_count_dtype = np.dtype([('alu_offset', np.uint16), ('count', np.uint32)])
    
    def __init__(self, num_16, data_obj, offset=None) :
        self.data_obj = data_obj
        self.offset = offset
        self.num_data = data_obj.num_data_from_num(num_16, self.offset)
        
    def offset_counts(self) :
        offsets, counts = np.unique(self.num_data['alu_offset'], return_counts=True)
        offset_counts = np.zeros(offsets.size, dtype=self.offset_count_dtype)
        offset_counts['alu_offset'] = offsets
        offset_counts['count'] = counts
        return offset_counts
    
    def num_repeat_indexes(self) :
        return np.unique(self.num_data['repeat_index'])
    
    def num_repeats_num_data(self) :
        ri = self.num_repeat_indexes()
        num_data = self.data_obj.num_data_from_repeat_indexes(ri)
        return num_data
        
    def num_repeats_offset_num_counts(self) :
        repeat_num_data = self.num_repeats_num_data()
        oao = offset_anal_cls(repeat_num_data, self.data_obj)
        oao.do_work()
        return oao.offset_num_counts
    
class num_assoc_anal_cls(object) :
    def __init__(self, num_1, offset_1, num_2, offset_2, data_obj) :
        self.num_1 = num_1
        self.offset_1 = offset_1
        self.num_2 = num_2
        self.offset_2 = offset_2
        self.data_obj = data_obj
        self.num_data = None
        
    def read_assoc_num_data(self) :
        ri_1 = self.data_obj.repeat_indexes_from_num(self.num_1, offset=self.offset_1)
        ri_2 = self.data_obj.repeat_indexes_from_num(self.num_2, offset=self.offset_2)
        ri = np.intersect1d(ri_1, ri_2)
        self.num_data = self.data_obj.num_data_from_repeat_indexes(ri)
        
    def assoc_offset_num_counts(self) :
        if self.num_data is None :
            self.read_assoc_num_data()
        oao = offset_anal_cls(self.num_data, self.data_obj)
        oao.do_work()
        return oao.offset_num_counts
        
class offset_pair_associations(object) :
    data_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('assoc_count', np.uint32), ('rel_expr', np.float32)])
    
    def __init__(self, target_offset, pair_num_16, pair_offset, data_obj, min_count=100) :
        self.target_offset = target_offset
        self.pair_num_16 = pair_num_16
        self.pair_offset = pair_offset
        self.data_obj = data_obj
        self.min_count = min_count
        self.target_num_data = self.data_obj.num_data_from_offset(self.target_offset)
        self.target_count = self.target_num_data.size
        self.pair_repeat_indexes = self.data_obj.repeat_indexes_from_num(self.pair_num_16, offset=self.pair_offset)
        self.pair_count = self.pair_repeat_indexes.size
        self.pred_factor = float(self.pair_count)/(self.target_count)
        
    def process_target(self) :
        self.target_num_data.sort(order=['num_16', 'repeat_index'])
        nums, starts, counts = np.unique(self.target_num_data['num_16'], return_index=True, return_counts=True)
        m = counts >= self.min_count
        nums = nums[m]
        starts = starts[m]
        counts = counts[m]
        out_data = []
        for i in range(nums.size) :
            num_16 = nums[i]
            start = starts[i]
            count = counts[i]
            bound = start + count
            num_repeat_indexes = self.target_num_data['repeat_index'][start:bound]
            assoc_repeat_indexes = np.intersect1d(num_repeat_indexes, self.pair_repeat_indexes)
            assoc_count = assoc_repeat_indexes.size
            pred_assoc = float(count)*self.pred_factor
            rel_expr = float(assoc_count)/pred_assoc
            out_data.append((num_16, count, assoc_count, rel_expr))
        self.association_data = np.array(out_data, self.data_dtype)


class standard_offset_associations(object) :
    standard_num_offset = ((26, 3277099518, 'c_26', 'r_26'), (117, 3338678528, 'c_117', 'r_117'),
                   (159, 3962915271, 'c_159', 'r_159'), (178, 2323162274, 'c_178', 'r_178'),  
                   (216, 2331586085, 'c_216', 'r_216'), (245, 3832882666, 'c_245', 'r_245'), 
                   (279, 3707764736, 'c_279', 'r_279'), (281, 0, 'c_280', 'r_280') )




                         
class offset_nums_anal_cls(object) :
    count_dtype = np.dtype([('num_16', np.uint32), ('start', np.uint32), ('count', np.uint32)])
    
    def __init__(self, offset, data_obj) :
        self.offset = offset
        self.data_obj = data_obj
        self.num_data = data_obj.num_data_from_offset(self.offset)
        self.num_counts = None

    def build_num_counts(self)  :
        self.num_data.sort(order=['num_16', 'repeat_index', 'alu_offset'])
        nums, starts, counts = np.unique(self.num_data['num_16'], return_index=True, return_counts=True)
        nc = np.zeros(nums.size, dtype=self.count_dtype)
        nc['num_16'] = nums
        nc['start'] = starts
        nc['count'] = counts
        self.num_counts = nc
        
    def gen_count_num_data(self, min_num_count=100) :
        if self.num_counts is None :
            self.build_num_counts()
        m = self.num_counts['count'] >= min_num_count
        nc = self.num_counts[m]
        for num_16, start, count in nc :
            bound = start + count
            nc_data = self.num_data[start:bound]
            yield num_16, count, nc_data


        