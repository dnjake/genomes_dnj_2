# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb
import genomes_dnj_2.autosome_snp_data.chrom_snp_data_rdr as sdr


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

class alu_snp_data_base_cls(object) :
    alu_sequence_data_folder = 'grch37_hg19_alu_data'
    data_folder = os.path.join(mod_dir, alu_sequence_data_folder)
    
    def __init__(self) :
        self.read_data()

    def read_data(self) :
        file_path = os.path.join(self.data_folder, self.file_name)
        h5 = tb.open_file(file_path, 'r')
        snp_data_table = getattr(h5.root, self.snp_data_table_name)
        self.snp_data = snp_data_table[:]
        alu_repeat_table = getattr(h5.root, self.alu_repeat_table_name)
        self.repeats = alu_repeat_table[:]
        h5.close()
            
    
    
class pos_alu_snp_data_cls(alu_snp_data_base_cls) :
    file_name = 'pos_alu_repeat_snps.h5'
    snp_data_table_name = 'pos_alu_repeat_snps'
    alu_repeat_table_name = 'alu_plus_strand_repeats'

    def __init__(self) :
        alu_snp_data_base_cls.__init__(self)

class neg_alu_snp_data_cls(alu_snp_data_base_cls) :
    file_name = 'neg_alu_repeat_snps.h5'
    snp_data_table_name = 'neg_alu_repeat_snps'
    alu_repeat_table_name = 'alu_minus_strand_repeats'
    
    def __init__(self) :
        alu_snp_data_base_cls.__init__(self)
    
    
class snp_rdr_cls(object) :
    data_dtype = np.dtype([('chrom', np.uint16), ('snp_index', np.uint32), ('snp_data', 'O')])
    
    
    def read_chrom_snps(self, chrom, snp_indexes) :
        rdr = sdr.chrom_snp_data_tables_cls(chrom)
        snps = rdr.snp_data_table[snp_indexes]
        rdr.close()
        return snps
        
    def read_snps(self, chroms, snp_indexes) :
        '''
        od = np.zeros(chroms.size, dtype=self.data_dtype)
        od['chrom'] = chroms
        od['snp_index'] = snp_indexes
        '''
        out_data = []
        u_chroms = np.unique(chroms)
        for chrom in u_chroms :
            m = chroms == chrom
            csd = self.read_chrom_snps(chrom, snp_indexes[m])
            out_data.append(csd)
        out_data = np.concatenate(out_data)
        return out_data
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        