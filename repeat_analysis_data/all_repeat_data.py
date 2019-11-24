# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

class all_repeats_cls(object) :
    repeat_file_name = 'indexed_sorted_repeats.h5'
    repeat_table_name = 'indexed_sorted_repeats'
    repeat_folder = 'grch37_hg19_annotation_data'
    
    def __init__(self) :
        self.read_repeat_data()

    def read_repeat_data(self) :
        folder = self.repeat_folder
        file_name = self.repeat_file_name
        local_path = os.path.join(folder, file_name) 
        file_path = os.path.join(mod_dir, local_path)
        h5 = tb.open_file(file_path, 'r')
        repeat_table = getattr(h5.root, self.repeat_table_name)
        self.repeat_data = repeat_table[:]
        h5.close()


class alu_repeat_cls(object) :
    unique_data_descr = [('assoc_num', np.uint32), ('num_starts', np.uint32), ('num_counts', np.uint32)]
    def __init__(self, strand=None) :
        self.strand = strand
        self.read_repeats()
        
    def read_repeats(self) :
        all_repeats_obj = all_repeats_cls()
        repeat_data = all_repeats_obj.repeat_data
        m_alu = repeat_data['repeat_family'] == 'alu'
        repeat_data = repeat_data[m_alu]
        self.repeat_data = repeat_data
        if self.strand is not None :
            m = repeat_data['strand'] == self.strand
            self.repeat_data = repeat_data[m]
        
    def chrom_alu_data(self, chrom) :
        chrom_bound = chrom + 1
        rdc = self.repeat_data['chrom'] 
        di = rdc.searchsorted((chrom, chrom_bound))
        return self.repeat_data[di[0]:di[1]]
    
    def repeats_from_indexes(self, repeat_indexes) :
        inds_rd = self.repeat_data['index'].searchsorted(repeat_indexes)
        return self.repeat_data[inds_rd]
    
    def repeats_from_name(self, repeat_name) :
        m = self.repeat_data['repeat_name'] == repeat_name
        return self.repeat_data[m]
    
    
class pos_alu_repeat_cls(alu_repeat_cls) :    
    def __init__(self) :
        alu_repeat_cls.__init__(self, strand='+')
        
class neg_alu_repeat_cls(alu_repeat_cls) :    
    def __init__(self) :
        alu_repeat_cls.__init__(self, strand='-')

