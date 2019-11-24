# -*- coding: utf-8 -*-

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

    def chrom_repeat_snps(self, chrom) :
        bound = chrom + 1
        inds = self.snp_data['chrom'].searchsorted([chrom, bound])
        return self.snp_data[inds[0]:inds[1]]
    
class pos_alu_snp_data_cls(alu_snp_data_base_cls) :
    file_name = 'pos_alu_repeat_full_snps.h5'
    snp_data_table_name = 'pos_alu_repeat_full_snps'
    alu_repeat_table_name = 'alu_plus_strand_repeats'

    def __init__(self) :
        alu_snp_data_base_cls.__init__(self)

class neg_alu_snp_data_cls(alu_snp_data_base_cls) :
    file_name = 'neg_alu_repeat_full_snps.h5'
    snp_data_table_name = 'neg_alu_repeat_full_snps'
    alu_repeat_table_name = 'alu_minus_strand_repeats'
    
    def __init__(self) :
        alu_snp_data_base_cls.__init__(self)


class dna_rdr_cls(object) :
    letters_folder = 'grch37_hg19_chrom_dna_letters'
    file_name_start = 'grch37_chr'
    file_name_end = '_fa.h5'
    array_name_start = 'grch37_chr_'

    def __init__(self, chrom) :
        self.chrom = chrom
        self.chrom_str = str(chrom)
        self.read_file()
                
    def read_file(self) :
        file_name = self.file_name_start + self.chrom_str + self.file_name_end
        local_path = os.path.join(self.letters_folder, file_name)
        file_path = os.path.join(mod_dir, local_path)
        array_name = self.array_name_start + self.chrom_str
        h5 = tb.open_file(file_path, 'r')
        array = getattr(h5.root, array_name)
        self.dna = array[:]
        h5.close()
    
    
    
    
    
    
    
    
    
    
    
    