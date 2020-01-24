# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d
import genomes_dnj_2.repeat_analysis.parital_sequences_analysis as psa


class match_atom(object) :
    def __init__(self, base_offset, match_delta_offset, match_seqs, match_repeat_indexes, test_seq_count) :
        self.base_offset = base_offset,
        self.match_delta_offset = match_delta_offset
        self.match_seqs = match_seqs
        self.match_repeat_indexes = match_repeat_indexes
        self.test_seq_count = test_seq_count
        
class base_offset_match_cls(object) :
    def __init__(self, base_offset, match_atoms, match_repeat_indexes, contrib_repeat_indexes, test_repeat_count) :
        self.base_offset = base_offset
        self.match_atoms = match_atoms
        self.match_repeat_indexes = match_repeat_indexes
        self.contrib_repeat_indexes = contrib_repeat_indexes
        self.test_repeat_count = test_repeat_count

class match_pattern_cls(object) :
    def __init__(self, delta_offset, match_seqs) :
        self.delta_offset = delta_offset
        self.match_seqs = match_seqs
        


class offset_cluster_match_cls(object) :
    
    def __init__(self, base_offsets, patterns, osq_obj) :
        self.base_offsets = base_offsets
        self.patterns = patterns
        self.osq_obj = osq_obj
        self.base_offset_matches = []
        self.match_repeat_indexes = None
        self.test_seq_count = 0
        
    def do_base_offset_match(self, base_offset) :
        match_atoms = []
        for p in self.patterns :
            match_offset = base_offset + p.delta_offset
            match_obj = self.osq_obj.offset_seq_match_obj(match_offset)
            p_ris = match_obj.ris_from_seqs(p.match_seqs)
            test_seq_count = match_obj.test_seq_count()
            match_atoms.append(match_atom(base_offset, p.delta_offset, p.match_seqs, p_ris, test_seq_count))
        ma = match_atoms[0]
        match_repeat_indexes = ma.match_repeat_indexes
        test_seq_count = ma.test_seq_count
        for ma in match_atoms[1:] :
            match_repeat_indexes = np.union1d(match_repeat_indexes, ma.match_repeat_indexes)
            if ma.test_seq_count > test_seq_count :
                test_seq_count = ma.test_seq_count
        if self.match_repeat_indexes is None :
            self.match_repeat_indexes = match_repeat_indexes
            contrib_repeat_indexes = match_repeat_indexes
        else :
            contrib_repeat_indexes = np.setdiff1d(match_repeat_indexes, self.match_repeat_indexes)
            self.match_repeat_indexes = np.union1d(self.match_repeat_indexes, match_repeat_indexes)
        bmo = base_offset_match_cls(base_offset, self.patterns, match_repeat_indexes, contrib_repeat_indexes, test_seq_count)
        self.base_offset_matches.append(bmo)
        if test_seq_count > self.test_seq_count :
            self.test_seq_count = test_seq_count
        
    def do_base_offsets(self) :
        for offset in self.base_offsets :
            self.do_base_offset_match(offset)

    def total_test_ris(self) :
        test_ris = None
        for offset in self.base_offsets :
            offset_obj = self.osq_obj.offset_seq_match_obj(offset)
            offset_ris = offset_obj.tested_ris()
            if offset_ris is None :
                continue
            if test_ris is None :
                test_ris = offset_ris
            else :
                test_ris = np.union1d(test_ris, offset_ris)
        return test_ris

class offset_num_16_match_cls(object) :
    def __init__(self, base_offsets, osq_obj) :
        self.base_offsets = base_offsets
        self.osq_obj = osq_obj
        
    def ris_num_16(self, delta_offset, num_16) :
        match_ris = None
        for base_offset in self.base_offsets :
            match_offset = base_offset + delta_offset
            match_obj = self.osq_obj.offset_seq_match_obj(match_offset)
            offset_ris = match_obj.num_16_repeat_indexes(num_16)
            if offset_ris is None :
                continue
            if match_ris is None :
                match_ris = offset_ris
            else :
                match_ris = np.union1d(match_ris, offset_ris)
        return match_ris

            
class cluster_match_num_counts_cls(object) :
    num_count_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('letters', 'S16')])
    
    def __init__(self, base_offset_matches, num_data_obj) :
        self.base_offset_matches = base_offset_matches
        self.num_data_obj = num_data_obj
        
    def build_match_num_counts(self, base_offset_delta=0) :
        match_num_data = []
        for bmo in self.base_offset_matches :
            num_offset = bmo.base_offset + base_offset_delta
            ond = self.num_data_obj.num_data_from_offset_and_repeat_indexes(num_offset, bmo.contrib_repeat_indexes)
            match_num_data.append(ond)
        match_num_data = np.concatenate(match_num_data)
        nums, counts = np.unique(match_num_data['num_16'], return_counts=True)
        mnc = np.zeros(nums.size, dtype=self.num_count_dtype)
        mnc['num_16'] = nums
        mnc['count'] = counts
        mnc['letters'] = n16d.num_16_strs(mnc['num_16'])
        self.num_counts = mnc
        
    def test_repeat_count(self) :
        bmo = self.base_offset_matches[0]
        count = bmo.test_repeat_count
        for bmo in self.base_offset_matches[1:] :
            if bmo.test_repeat_count > count :
                count = bmo.test_repeat_count
        return count
        
    def build_sorted_match_num_counts(self, base_offset_delta=0) :
        self.build_match_num_counts(base_offset_delta)
        self.num_counts.sort(order='count')
        self.num_counts = self.num_counts[::-1]
                                          
                                          
class cluster_match_seq_counts_cls(object) :
    seq_count_dtype = np.dtype([('letters', 'S16'), ('count', np.uint32)])
    
    def __init__(self, base_offset_matches, num_data_obj) :
        self.base_offset_matches = base_offset_matches
        self.num_data_obj = num_data_obj

    def test_repeat_count(self) :
        bmo = self.base_offset_matches[0]
        count = bmo.test_repeat_count
        for bmo in self.base_offset_matches[1:] :
            if bmo.test_repeat_count > count :
                count = bmo.test_repeat_count
        return count

    def build_match_seq_counts(self, seq_size, base_offset_delta=0) :
        match_nums = []
        for bmo in self.base_offset_matches :
            num_offset = bmo.base_offset + base_offset_delta
            ond = self.num_data_obj.num_data_from_offset_and_repeat_indexes(num_offset, bmo.contrib_repeat_indexes)
            mn = psa.masked_nums(ond['num_16'], seq_size)
            match_nums.append(mn)
        match_nums = np.concatenate(match_nums)
        nums, counts = np.unique(match_nums, return_counts=True)
        mnc = np.zeros(nums.size, dtype=self.seq_count_dtype)
        mnc['letters'] = n16d.num_16_short_strs(nums, seq_size)
        mnc['count'] = counts
        self.seq_counts = mnc
        
    def build_sorted_match_seq_counts(self, seq_size, base_offset_delta=0) :
        self.build_match_seq_counts(seq_size, base_offset_delta)
        self.seq_counts.sort(order='count')
        self.seq_counts = self.seq_counts[::-1]      
                                          
class cluster_seq_match_ris_cls(object) :
    
    def __init__(self, base_offset_matches, num_data_obj) :
        self.base_offset_matches = base_offset_matches
        self.num_data_obj = num_data_obj

    def test_repeat_count(self) :
        bmo = self.base_offset_matches[0]
        count = bmo.test_repeat_count
        for bmo in self.base_offset_matches[1:] :
            if bmo.test_repeat_count > count :
                count = bmo.test_repeat_count
        return count

    def seq_match_ris(self, match_seq, base_offset_delta=0) :
        match_ris = []
        seq_size = len(match_seq)
        match_seq_num = psa.num_16_from_substr(match_seq)
        for bmo in self.base_offset_matches :
            num_offset = bmo.base_offset + base_offset_delta
            ond = self.num_data_obj.num_data_from_offset_and_repeat_indexes(num_offset, bmo.contrib_repeat_indexes)
            mn = psa.masked_nums(ond['num_16'], seq_size)
            m_match = mn == match_seq_num
            if m_match.sum() > 0 :
                offset_ris = ond['repeat_index'][m_match]
                match_ris.append(offset_ris)
        if len(match_ris) > 0 :
            match_ris = np.concatenate(match_ris)
            return match_ris
        
    def or_seq_match_ris(self, offset_seqs) :
        base_offset_delta, match_seqs = offset_seqs
        match_ris = None
        for seq in match_seqs :
            ris = self.seq_match_ris(seq, base_offset_delta)
            if ris is None :
                continue
            if match_ris is None :
                match_ris = ris
            else :
                match_ris = np.union1d(match_ris, ris)
        return match_ris
    
    def or_pattern_match_ris(self, patterns) :
        match_ris = None
        for offset_seqs in patterns :
            ris = self.or_seq_match_ris(offset_seqs)
            if ris is None :
                continue
            if match_ris is None :
                match_ris = ris
            else :
                match_ris = np.union1d(match_ris, ris)
        return match_ris
    
    def and_pattern_match_ris(self, patterns) :
        match_ris = None
        for offset_seqs in patterns :
            ris = self.or_seq_match_ris(offset_seqs)
            if match_ris is None :
                match_ris = ris
            else :
                match_ris = np.intersect1d(match_ris, ris)
        return match_ris
            
         
        
class no_match_num_counts(object) :
    num_count_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('letters', 'S16')])
    def __init__(self, match_ris, num_data_obj) :
        self.match_ris = match_ris
        self.num_data_obj = num_data_obj
        all_ris = num_data_obj.data_all_repeat_indexes()
        self.no_match_ris = np.setdiff1d(all_ris, self.match_ris)
        
    def offset_num_counts(self, offset) :
        no_match_num_data = self.num_data_obj.num_data_from_offset_and_repeat_indexes(offset, self.no_match_ris)
        nums, counts = np.unique(no_match_num_data['num_16'], return_counts=True)
        mnc = np.zeros(nums.size, dtype=self.num_count_dtype)
        mnc['num_16'] = nums
        mnc['count'] = counts
        mnc['letters'] = n16d.num_16_strs(mnc['num_16'])
        return mnc

    def sorted_offset_num_counts(self, offset) :
        num_counts = self.offset_num_counts(offset)
        num_counts.sort(order='count')
        num_counts = num_counts[::-1]
        return num_counts
        
int_fmt =  '{:d}'
big_fmt =  '{:,}'
float_fmt = '{:.2f}'
exp_fmt = '{:.1e}'
            

class alignment_run_cls(object) :

    def __init__(self, align_obj, osqo, ndo) :
        self.align_obj = align_obj
        self.segment = align_obj.segment
        self.alu_offset = align_obj.alu_offset
        self.segment_offset = align_obj.segment_offset
        self.osqo = osqo
        self.ndo = ndo
        self.ocmo = offset_cluster_match_cls(align_obj.base_offsets, align_obj.patterns, osqo)
        self.ocmo.do_base_offsets()
        self.aligned_ris = self.ocmo.match_repeat_indexes
        self.test_seq_count = self.ocmo.test_seq_count
        frac_aligned = float(self.aligned_ris.size)/float(self.test_seq_count)
        print('alu offset:', align_obj.alu_offset, '\nsegment:', self.align_obj.segment,
              '\nsegment offset:', align_obj.segment_offset)
        print('test sequence count', big_fmt.format(self.test_seq_count))
        print('aligned repeat count', big_fmt.format(self.aligned_ris.size))
        print('fraction aligned', float_fmt.format(frac_aligned))
        
    def get_top_num_16(self, base_offset_delta, num_count=100) :
        cmnco = cluster_match_num_counts_cls(self.ocmo.base_offset_matches, self.ndo)
        cmnco.build_sorted_match_num_counts(base_offset_delta)
        nc = cmnco.num_counts
        return nc[:num_count]
        
    def top_100_num_16(self, base_offset_delta, num_count=100) :
        cmnco = cluster_match_num_counts_cls(self.ocmo.base_offset_matches, self.ndo)
        cmnco.build_sorted_match_num_counts(base_offset_delta)
        nc = cmnco.num_counts
        total_aligned = nc['count'].sum()
        top_100 = nc[:num_count]
        top_100_sum = top_100['count'].sum()
        top_100_frac = float(top_100_sum)/float(total_aligned)
        alu_offset = self.alu_offset + base_offset_delta
        segment_offset = self.segment_offset + base_offset_delta
        segment = self.segment
        print('alu offset:', alu_offset, '\nsegment:', segment, '\nsegment offset:', segment_offset)
        base_offsets = self.align_obj.base_offsets
        print('base offsets', base_offsets, '\nsequence count', big_fmt.format(nc.size),
              '\naligned repeat count', big_fmt.format(total_aligned))
        num_count_str = 'top ' + str(num_count)
        repeat_str = num_count_str + ' repeat_count'
        frac_str = '\n' + num_count_str + ' fraction'
        base_count_str = num_count_str + ' 16 base_counts\n'
        print(repeat_str, big_fmt.format(top_100_sum), frac_str, float_fmt.format(top_100_frac))
        count_offsets = [offset+base_offset_delta for offset in base_offsets]
        print('count offsets', count_offsets)
        print(base_count_str, top_100)

    def top_100_seqs(self, seq_offset, seq_size) :
        smnco = cluster_match_seq_counts_cls(self.ocmo.base_offset_matches, self.ndo)
        smnco.build_sorted_match_seq_counts(seq_size, seq_offset)
        sc = smnco.seq_counts
        alu_offset = self.alu_offset + seq_offset
        segment_offset = self.segment_offset + seq_offset
        segment = self.segment
        print('alu offset:', alu_offset, '\nsegment:', segment, '\nsegment offset:', segment_offset)
        base_offsets = self.align_obj.base_offsets
        print('base offsets', base_offsets, '\nsequence count', big_fmt.format(sc.size),
              '\nrepeat count', big_fmt.format(sc['count'].sum()))
        seq_offsets = [offset+seq_offset for offset in base_offsets]
        print('sequence offsets', seq_offsets)
        print('top 100 sequence counts\n', sc[:100])








        
        
