# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d


dna_base_codes = (('A', 0), ('C', 1), ('G', 2), ('T', 3))
a_str_a = np.empty(16, dtype='S1')
a_str_a[:] = 'A'
char_mults = np.array([4**15, 4**14, 4**13, 4**12, 4**11, 4**10, 4**9, 4**8, 4**7, 4**6,
                       4**5, 4**4, 4**3, 4**2, 4**1, 4**0], dtype=np.uint32)
    


def masked_nums(nums, mask_size) :
    index = mask_size - 1
    mask_mult = char_mults[index]
    masked_nums = np.floor_divide(nums, mask_mult)
    masked_nums*= mask_mult
    return masked_nums

def num_16_from_substr(letters) :
    num_a = a_str_a.copy()
    letter_a = np.frombuffer(letters, dtype='S1')
    s = letter_a.size
    num_a[:s] = letter_a
    out_a = np.empty(num_a.size, dtype=np.uint32)
    for c_l, c_n in dna_base_codes :
        m = num_a == c_l
        out_a[m] = c_n
    out_a*=char_mults
    return out_a.sum()

masked_cg = num_16_from_substr('CG')
masked_tg = num_16_from_substr('TG')
masked_ca = num_16_from_substr('CA')
masked_ttagc = num_16_from_substr('TTAGC')
tail_seqs = ['CGTCTC', 'TGTCTC', 'CATCTC']


class masked_num_match_cls(object) :
    def __init__(self, num_16s, mask_size) :
        self.mask_size = mask_size
        self.masked_nums = masked_nums(num_16s, mask_size)
        
    def match_num(self, masked_num_16) :
        return self.masked_nums == masked_num_16

    
class cg_and_cg_derived_counts_cls(object) :
    offset_count_dtype = np.dtype([('offset', np.uint16), ('count', np.uint32)])
    @classmethod
    def ris_with_cgd_match(cls, offset, num_data) :
        m = num_data['alu_offset'] == offset
        num_data = num_data[m]
        masked_num_match_obj = masked_num_match_cls(num_data['num_16'], 2)
        m = masked_num_match_obj.match_num(masked_cg)
        ris = num_data['repeat_index'][m]
        m = masked_num_match_obj.match_num(masked_tg)
        ris_m = num_data['repeat_index'][m]
        ris = np.union1d(ris, ris_m)
        m = masked_num_match_obj.match_num(masked_ca)
        ris_m = num_data['repeat_index'][m]
        ris = np.union1d(ris, ris_m)
        return ris

    @classmethod
    def from_repeat_indexes(cls, repeat_indexes, data_obj) :
        #num_data = data_obj.num_data_from_repeat_indexes(repeat_indexes)
        return cls(repeat_indexes, data_obj)
    
    @classmethod
    def from_offsets(cls, offsets, data_obj):
        num_data = data_obj.num_data
        ris = None
        for offset in offsets :
            offset_ris = cls.ris_with_cgd_match(offset, num_data)
            if ris is None :
                ris = offset_ris
            else :
                ris = np.intersect1d(ris, offset_ris)
        return cls(ris, data_obj)
    
    @classmethod
    def from_yes_no_offsets(cls, data_obj, yes_offsets=None, no_offsets=None) :
        ris = None
        num_data = data_obj.num_data
        if yes_offsets is None :
            ris = data_obj.data_all_repeat_indexes()
        else :
            for offset in yes_offsets :
                offset_ris = cls.ris_with_cgd_match(offset, num_data)
                if ris is None :
                    ris = offset_ris
                else :
                    ris = np.intersect1d(ris, offset_ris)
        if no_offsets is not None :
            no_ris = None
            for offset in no_offsets :
                offset_ris = cls.ris_with_cgd_match(offset, num_data)
                if no_ris is None :
                    no_ris = offset_ris
                else :
                    no_ris = np.union1d(no_ris, offset_ris)
            ris = np.setdiff1d(ris, no_ris)
        return cls(ris, data_obj)
        
    def __init__(self, repeat_indexes, data_obj) :
        self.repeat_indexes = repeat_indexes
        self.data_obj = data_obj
        self.repeats_num_data = data_obj.num_data_from_repeat_indexes(repeat_indexes)
        self.masked_num_match_obj = masked_num_match_cls(self.repeats_num_data['num_16'], 2)
        self.do_matches()

    def do_unique_counts(self, match_offsets) :
        offsets, counts = np.unique(match_offsets, return_counts=True)
        oc = np.zeros(offsets.size, self.offset_count_dtype)
        oc['offset'] = offsets
        oc['count'] = counts
        return oc

    def do_matches(self) :
        m_match = self.masked_num_match_obj.match_num(masked_cg)
        cg_match_offsets = self.repeats_num_data['alu_offset'][m_match]
        self.cg_offset_counts = self.do_unique_counts(cg_match_offsets)
        m_match = self.masked_num_match_obj.match_num(masked_tg)
        tg_match_offsets = self.repeats_num_data['alu_offset'][m_match]
        m_match = self.masked_num_match_obj.match_num(masked_ca)
        ca_match_offsets = self.repeats_num_data['alu_offset'][m_match]
        cg_derived_offsets = np.concatenate([cg_match_offsets, tg_match_offsets, ca_match_offsets])
        self.cg_derived_offset_counts = self.do_unique_counts(cg_derived_offsets)
        



class short_sequence_associations_cls(object) :
    rel_num_count_dtype = np.dtype([('num_16', np.uint32), ('repeat_offset_num_count', np.uint32), 
                                    ('data_offset_num_count', np.uint32), ('rel_expr', np.float32), ('letters', 'S16')])
    

    @staticmethod
    def repeat_indexes_from_offset_seq(offset, seq, data_obj) :
        num_data = data_obj.num_data
        m = num_data['alu_offset'] == offset
        seq_num = num_16_from_substr(seq)
        offset_nums = masked_nums(num_data['num_16'][m], len(seq))
        m_nums = offset_nums == seq_num
        repeat_indexes = num_data['repeat_index'][m]
        repeat_indexes = repeat_indexes[m_nums]
        return repeat_indexes
    
    @classmethod
    def from_offset_seq_pairs(cls, offset_seqs, data_obj) :
        repeat_indexes = None
        for offset, seq in offset_seqs :
            offset_repeat_indexes = cls.repeat_indexes_from_offset_seq(offset, seq, data_obj)
            if repeat_indexes is None :
                repeat_indexes = offset_repeat_indexes
            else :
                repeat_indexes = np.intersect1d(repeat_indexes, offset_repeat_indexes)
        return cls(repeat_indexes, data_obj)
    

    @classmethod
    def from_offset_seq(cls, offset, seq, data_obj) :
        num_data = data_obj.num_data
        m = num_data['alu_offset'] == offset
        seq_num = num_16_from_substr(seq)
        offset_nums = masked_nums(num_data['num_16'][m], len(seq))
        m_nums = offset_nums == seq_num
        repeat_indexes = num_data['repeat_index'][m]
        repeat_indexes = repeat_indexes[m_nums]
        return cls(repeat_indexes, data_obj)
                
    def __init__(self, repeat_indexes, data_obj) :
        self.repeat_indexes = repeat_indexes
        self.data_obj = data_obj
        self.total_num_data_obj = data_obj.num_data_obj
        self.repeats_num_data_obj = self.data_obj.obj_from_repeat_indexes(self.repeat_indexes)
        
    def counts_for_offset_and_seq_size(self, offset, seq_size, num_data) :
        m = num_data['alu_offset'] == offset
        seq_nums = masked_nums(num_data['num_16'][m], seq_size)
        nums, counts = np.unique(seq_nums, return_counts=True)
        return nums, counts
    
    def insert_num_letters(self, str_size, count_data) :
        count_data['letters'] = n16d.num_16_short_strs(count_data['num_16'], str_size)
        
    def merge_repeat_num_counts(self, nums, counts, count_data) :
        inds = count_data['num_16'].searchsorted(nums)
        m = inds < count_data.size
        m[m] = nums[m] == count_data['num_16'][inds[m]]
        count_data['repeat_offset_num_count'][inds[m]] = counts[m]
    
    def calc_rel_expr(self, repeat_count, total_count, count_data) :
        pred_factor = float(repeat_count)/float(total_count)
        pred_repeat_count = pred_factor*count_data['data_offset_num_count'].astype(np.float32)        
        count_data['rel_expr'] = count_data['repeat_offset_num_count'].astype(np.float32)/pred_repeat_count

    def rel_counts_for_offset_and_seq_size(self, offset, seq_size, min_count=100) :
        num_data = self.total_num_data_obj.num_data
        nums, counts = self.counts_for_offset_and_seq_size(offset, seq_size, num_data)
        total_offset_count = counts.sum()
        m = counts >= min_count
        data_size = m.sum()
        cd = np.zeros(data_size, dtype=self.rel_num_count_dtype)
        cd['num_16'] = nums[m]
        cd['data_offset_num_count'] = counts[m]
        num_data = self.repeats_num_data_obj.num_data
        nums, counts = self.counts_for_offset_and_seq_size(offset, seq_size, num_data)
        repeat_offset_count = counts.sum()
        self.merge_repeat_num_counts(nums, counts, cd)
        self.calc_rel_expr(repeat_offset_count, total_offset_count, cd)
        self.insert_num_letters(seq_size, cd)
        cd.sort(order='data_offset_num_count')
        cd = cd[::-1]
        return cd
    
    def repeat_sorted_rel_seq_counts(self, offset, seq_size, min_count=100) :
        rsc = self.rel_counts_for_offset_and_seq_size(offset, seq_size, min_count)
        rsc.sort(order=['repeat_offset_num_count', 'data_offset_num_count'])
        rsc = rsc[::-1]
        return rsc
        
class subseq_match_cls(object) :
    offset_count_dtype = np.dtype([('offset', np.uint16), ('count', np.uint32)])
    
    def __init__(self, num_data, mask_size) :
        self.num_data = num_data
        self.mask_size = mask_size
        self.seq_match_obj = masked_num_match_cls(num_data['num_16'], mask_size)
        
    def seq_offsets(self, seq) :
        assert len(seq) == self.mask_size
        masked_seq = num_16_from_substr(seq)
        m_match = self.seq_match_obj.match_num(masked_seq)
        return self.num_data['alu_offset'][m_match]
    
    def offset_seq_repeat_indexes(self, offset, seq) :
        assert len(seq) == self.mask_size
        masked_seq = num_16_from_substr(seq)
        m_match = self.seq_match_obj.match_num(masked_seq)
        m_offset = self.num_data['alu_offset'] == offset
        m_match = np.logical_and(m_match, m_offset)
        return self.num_data['repeat_index'][m_match]
        
    def offset_counts_for_seqs(self, seqs) :
        if type(seqs) == str :
            seqs = [seqs]
        offsets = []            
        for seq in seqs :
            offsets.append(self.seq_offsets(seq))
        offsets = np.concatenate(offsets)
        offsets, counts = np.unique(offsets, return_counts=True)
        od = np.zeros(offsets.size, dtype=self.offset_count_dtype)
        od['offset'] = offsets
        od['count'] = counts
        return od
    
    def repeat_indexes_for_offset_seqs(self, offset, seqs) :        
        if type(seqs) == str :
            seqs = [seqs]
        repeat_indexes = []            
        for seq in seqs :
            repeat_indexes.append(self.offset_seq_repeat_indexes(offset, seq))
        repeat_indexes = np.concatenate(repeat_indexes)
        repeat_indexes.sort()
        return repeat_indexes
    
    
    