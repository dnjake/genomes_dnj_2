# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np
import genomes_dnj_2.repeat_analysis.parital_sequences_analysis as psa
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d

'''
It appears I am repeatedly fetching num data for an offset to try a set of matches on
sequences at the offset.  I also want to know the total number of repeats that have
sequences to match at the offset.  I think I want a separate object to do the matching
for an offset.  That object will hold the offset num data and do tests on it.  A call
from the using code can get the object and do a number of matches on it.  It can also
get the test_repeat_count from the object
'''

seqs_tail = (('CGTCTC'), ('TGTCTC'), ('CATCTC'))
seqs_cg_derived = (('CG'), ('TG'), ('CA'))
seq_ttagc = 'TTAGC'

class seq_match_cls(object) :
    def __init__(self, num_data) :
        self.num_data = num_data

    def test_seq_count(self) :
        return self.num_data.size

    def sub_seq_repeat_indexes(self, seq) :
        seq_size = len(seq)
        smo = psa.masked_num_match_cls(self.num_data['num_16'], seq_size)
        seq_value = psa.num_16_from_substr(seq)
        m = smo.match_num(seq_value)
        ris = self.num_data['repeat_index'][m]
        return ris
    
    def ris_from_seqs(self, seqs) :
        ris = None
        for seq in seqs :
            seq_ri = self.sub_seq_repeat_indexes(seq)
            if ris is None :
                ris = seq_ri
            else :
                ris = np.union1d(ris, seq_ri)
        return ris
    

class offset_sequence_cls(object) :
    
    no_letters_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32)])
    num_count_dtype = np.dtype([('num_16', np.uint32), ('count', np.uint32), ('letters', 'S16')])
    
    @classmethod
    def from_strand_data_obj(cls, data_obj) :
        return cls(data_obj.num_data_obj)

    @classmethod
    def from_repeat_indexes(cls, repeat_indexes, data_obj) :
        num_data_obj = data_obj.obj_from_repeat_indexes(repeat_indexes)
        return cls(num_data_obj)

    @classmethod
    def from_offset_num_16(cls, offset, num_16, data_obj) :
        ndo = data_obj.num_data_obj
        nd_offset = ndo.num_data_from_num(num_16, offset)
        ris_offset = nd_offset['repeat_index']
        return  cls.from_repeat_indexes(ris_offset, data_obj)

    @classmethod
    def from_num_16(cls, num_16, data_obj) :
        ndo = data_obj.num_data_obj
        nd_num = ndo.num_data_from_num(num_16)
        ris_num = nd_num['repeat_index']
        return  cls.from_repeat_indexes(ris_num, data_obj)

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
        return cls.from_repeat_indexes(repeat_indexes, data_obj) 

    @classmethod
    def from_alu_name(cls, alu_name, data_obj) :
        m = data_obj.repeats['repeat_name'] == alu_name
        repeat_indexes = data_obj.repeats['index'][m]
        return cls.from_repeat_indexes(repeat_indexes, data_obj)
    
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
        return cls.from_repeat_indexes(repeat_indexes, data_obj)
    
    def __init__(self, num_data_obj) :
        self.num_data_obj = num_data_obj

    def offset_sub_seq_repeat_indexes(self, offset, seq) :
        offset_num_data = self.num_data_obj.num_data_from_offset(offset)
        seq_size = len(seq)
        smo = psa.masked_num_match_cls(offset_num_data['num_16'], seq_size)
        seq_value = psa.num_16_from_substr(seq)
        m = smo.match_num(seq_value)
        ris = offset_num_data['repeat_index'][m]
        return ris

    def ris_from_offset_seqs(self, offset, seqs) :
        ris = None
        for seq in seqs :
            seq_ri = self.offset_sub_seq_repeat_indexes(offset, seq)
            if ris is None :
                ris = seq_ri
            else :
                ris = np.union1d(ris, seq_ri)
        return ris
    
    def offset_seq_match_obj(self, offset) :
        offset_num_data = self.num_data_obj.num_data_from_offset(offset)
        return seq_match_cls(offset_num_data)

    
    def ris_from_and_offset_seqs_pairs(self, pairs) :
        ris = None
        for offset, seqs in pairs :
            offset_ris = self.ris_from_offset_seqs(offset, seqs)
            if  (ris is None) or (offset_ris is None) :
                ris = offset_ris
            else :
                ris = np.intersect1d(ris, offset_ris)
        return ris

    def offset_num_16_repeat_indexes(self, offset, num_16) :
        ndo = self.num_data_obj
        ris = ndo.repeat_indexes_from_num(num_16, offset=offset)
        return ris
    
    def num_data_from_offset_and_repeat_indexes(self, offset, repeat_indexes) :
        ndo = self.num_data_obj
        return ndo.num_data_from_offset_and_repeat_indexes(offset, repeat_indexes)
    
    def num_counts_from_num_data(self, num_data, with_letters=False) :
        nums, counts = np.unique(num_data['num_16'], return_counts=True)
        if with_letters :
            nc_dtype = self.num_count_dtype
        else :
            nc_dtype = self.no_letters_dtype
        ncd = np.zeros(nums.size, dtype=nc_dtype)
        ncd['num_16'] = nums
        ncd['count'] = counts
        if with_letters :
            ncd['letters'] = n16d.num_16_strs(ncd['num_16'])
        return ncd
    
    def sorted_num_counts_from_num_data(self, num_data) :
        nc = self.num_counts_from_num_data(num_data, with_letters=True)
        nc.sort(order='count')
        nc = nc[::-1]
        return nc

    def sorted_num_counts_from_offset(self, offset) :
        nd = self.num_data_obj.num_data_from_offset(offset)
        return self.sorted_num_counts_from_num_data(nd)

    def seq_counts_from_num_data(self, num_data, seq_size, with_letters=False) :
        masked_nums = psa.masked_nums(num_data['num_16'], seq_size)
        nums, counts = np.unique(masked_nums, return_counts=True)
        if with_letters :
            nc_dtype = self.num_count_dtype
        else :
            nc_dtype = self.no_letters_dtype
        ncd = np.zeros(nums.size, dtype=nc_dtype)
        ncd['num_16'] = nums
        ncd['count'] = counts
        if with_letters :
            ncd['letters'] = n16d.num_16_short_strs(ncd['num_16'], seq_size)
        return ncd

    def sorted_seq_counts_from_num_data(self, num_data, seq_size) :
        nc = self.seq_counts_from_num_data(num_data, seq_size, with_letters=True)
        nc.sort(order='count')
        nc = nc[::-1]
        return nc

    def sorted_seq_counts_from_offset(self, offset, seq_size) :
        nd = self.num_data_obj.num_data_from_offset(offset)
        return self.sorted_seq_counts_from_num_data(nd, seq_size)

    def offset_value_matches(self, offset, match_nums, match_repeat_indexes=None) :
        ndo = self.num_data_obj
        if match_repeat_indexes is None :
            offset_num_data = ndo.num_data_from_offset(offset)
        else :
            offset_num_data = ndo.num_data_from_offset_and_repeat_indexes(offset, match_repeat_indexes)
        inds_matches = match_nums.searchsorted(offset_num_data['num_16'])
        m = inds_matches < match_nums.size
        inds_matches = inds_matches[m]
        inds_data = np.where(m)[0]
        m = match_nums[inds_matches] == offset_num_data['num_16'][inds_data]
        inds_data = inds_data[m]
        return offset_num_data[inds_data]

    def offset_sub_seq_matches(self, offset, match_nums, seq_size, match_repeat_indexes=None) :
        ndo = self.num_data_obj
        if match_repeat_indexes is None :
            offset_num_data = ndo.num_data_from_offset(offset)
        else :
            offset_num_data = ndo.num_data_from_offset_and_repeat_indexes(offset, match_repeat_indexes)
        offset_masked_nums = psa.masked_nums(offset_num_data['num_16'], seq_size)
        match_masked_nums = psa.masked_nums(match_nums, seq_size)
        match_masked_nums = np.unique(match_masked_nums)
        inds_matches = match_masked_nums.searchsorted(offset_masked_nums)
        m = inds_matches < match_masked_nums.size
        inds_matches = inds_matches[m]
        inds_data = np.where(m)[0]
        m = match_masked_nums[inds_matches] == offset_masked_nums[inds_data]
        inds_data = inds_data[m]
        return offset_num_data[inds_data]
    
    def most_common_num_16_data(self, offset, repeat_indexes=None) :
        ndo = self.num_data_obj
        if repeat_indexes is None :
            num_data = ndo.num_data_from_offset(offset)
        else :
            num_data = ndo.num_data_from_offset_and_repeat_indexes(offset, repeat_indexes)
        ncd = self.num_counts_from_num_data(num_data)
        ind_max = ncd['count'].argmax()
        max_num, max_count = ncd[ind_max]
        total_count = ncd['count'].sum()
        letters = n16d.str_from_num_16(max_num)
        return offset, letters, max_count, total_count, max_num
        
    def most_common_num_16_value(self, offset, repeat_indexes=None) :
        mcd = self.most_common_num_16_data(offset, repeat_indexes)
        return mcd[4]

    def most_common_seq(self, offset, seq_size, repeat_indexes=None) :
        ndo = self.num_data_obj
        if repeat_indexes is None :
            num_data = ndo.num_data_from_offset(offset)
        else :
            num_data = ndo.num_data_from_offset_and_repeat_indexes(offset, repeat_indexes)
        ncd = self.seq_counts_from_num_data(num_data, seq_size)
        ind_max = ncd['count'].argmax()
        max_num, max_count = ncd[ind_max]
        total_count = ncd['count'].sum()
        letters = n16d.str_from_num_16(max_num)
        letters = letters[:seq_size]
        return offset, letters, max_count, total_count, max_num
    
    '''
    def num_data_from_offset_and_nums(self, offset, match_num, seq_size=None) :
        offset_num_data = self.num_data_obj.num_data_from_offset(offset)
        if seq_size is None :
            data_num_16 = offset_num_data['num_16']
        else :
            data_num_16 = psa.masked_nums(offset_num_data['num_16'], seq_size)
        inds_match = match_num.searchsorted(data_num_16)
        m = inds_match < match_num.size
        inds_match = inds_match[m]
        inds_data = np.where(m)[0]
        m = match_num[inds_match] == data_num_16[inds_data]
        inds_data = inds_data[m]
        return offset_num_data[inds_data]
    '''

    def full_num_counts_from_ris_subset(self, offset, repeat_indexes, min_count=10, with_letters=False) :
        ris_num_data = self.num_data_from_offset_and_repeat_indexes(offset, repeat_indexes)
        ris_num_counts = self.num_counts_from_num_data(ris_num_data)
        m = ris_num_counts['count'] >= min_count
        ris_num_counts = ris_num_counts[m]
        full_num_data = self.offset_value_matches(offset, ris_num_counts['num_16'])
        full_num_counts = self.num_counts_from_num_data(full_num_data, with_letters)
        return full_num_counts
        
    def sorted_full_num_counts_from_ris_subset(self, offset, repeat_indexes, min_count=10) :
        full_num_counts = self.full_num_counts_from_ris_subset(offset, repeat_indexes, min_count, True)
        full_num_counts.sort(order='count')
        full_num_counts = full_num_counts[::-1]
        return full_num_counts
    
    def full_seq_counts_from_ris_subset(self, offset, repeat_indexes, seq_size, min_count=10, with_letters=False) :
        ris_num_data = self.num_data_from_offset_and_repeat_indexes(offset, repeat_indexes)
        ris_seq_counts = self.seq_counts_from_num_data(ris_num_data, seq_size)
        m = ris_seq_counts['count'] >= min_count
        ris_seq_counts = ris_seq_counts[m]
        full_num_data = self.offset_sub_seq_matches(offset, ris_seq_counts['num_16'], seq_size)
        full_seq_counts = self.seq_counts_from_num_data(full_num_data, seq_size, with_letters)
        return full_seq_counts
    
    def sorted_full_seq_counts_from_ris_subset(self, offset, repeat_indexes, seq_size, min_count=10) :
        full_seq_counts = self.full_seq_counts_from_ris_subset(offset, repeat_indexes, seq_size, min_count, True)
        full_seq_counts.sort(order='count')
        full_seq_counts = full_seq_counts[::-1]
        return full_seq_counts


num_tag = '<td style="text-align: right;">'
center_tag = '<td style="text-align: center;">'
header_tag = '<th style="text-align:center;">'
int_fmt =  '{:d}'
big_fmt =  '{:,}'

class offset_count_html_cls(object) :
    count_table_dtype = np.dtype([('offset', np.uint16), ('count', np.uint32), ('row', np.uint16), ('col', np.uint16)])
    row_count_dtype = np.dtype([('row', np.uint32), ('start', np.uint32), ('count', np.uint32)])
    display_size = 300
    row_size = 10
    cols = np.arange(row_size, dtype=np.uint32)
    
    def __init__(self, offset_counts) :
        self.offset_counts = offset_counts

    def build_count_table(self) :
        display_counts = np.zeros(self.display_size, dtype=self.count_table_dtype)
        display_counts['offset'] = np.arange(self.display_size)
        inds_display = display_counts['offset'].searchsorted(self.offset_counts['offset'])
        m = inds_display < display_counts.size
        inds_display = inds_display[m]
        inds_counts = np.where(m)[0]
        display_counts['count'][inds_display] = self.offset_counts['count'][inds_counts]
        display_counts['row'] = display_counts['offset']//self.row_size
        display_counts['col'] = display_counts['offset']%self.row_size
        m = display_counts['count'] > 0
        inds = np.where(m)[0]
        ind_max = inds.max()
        display_counts = display_counts[:ind_max]
        self.offset_counts = display_counts
        rows, starts, counts = np.unique(self.offset_counts['row'], return_index=True, return_counts=True)
        self.row_index = np.array(zip(rows,starts,counts), dtype=self.row_count_dtype)

    def init_table_html(self) :
        self.table_html = ['<table>\n']
        
    def build_header_html(self) :
        header_html = ['<thead>\n', '<tr>']
        header_html.append(header_tag)
        header_html.append('offset')
        header_html.append('</th>')
        for col in self.cols :
            header_html.append(header_tag)
            header_html.append(int_fmt.format(col))
            header_html.append('</th>')            
        header_html.append('</tr>\n</thead>\n')
        header_html = ''.join(header_html)
        self.table_html.append(header_html)
        
    def build_row_html(self, row, row_offset_counts) :
        row_html = ['<tr>']
        row_html.append(num_tag)
        row_val = row*self.row_size
        row_html.append(int_fmt.format(row_val))
        row_html.append('</td>')
        if row == 0 :
            row_html.append('<td></td>')
            row_offset_counts = row_offset_counts[1:]
        counts = row_offset_counts['count']
        for count in counts :
            row_html.append(num_tag)
            row_html.append(big_fmt.format(count))
            row_html.append('</td>')
        row_html.append('</tr>\n')
        row_html = ''.join(row_html)
        return row_html
        
                
    def build_rows(self) :
        for row, start, count in self.row_index :
            bound = start + count
            row_offset_counts = self.offset_counts[start:bound]
            row_html = self.build_row_html(row, row_offset_counts)
            self.table_html.append(row_html)
        
    def complete_table(self) :
        self.table_html.append('</table>\n')
        self.table_html = ''.join(self.table_html)

    def do_table(self) :
        self.build_count_table()
        self.init_table_html()
        self.build_header_html()
        self.build_rows()
        self.complete_table()


