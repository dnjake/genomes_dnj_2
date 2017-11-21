# -*- coding: utf-8 -*-
import tables as tb
import numpy as np
import os


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)


class chrom_snp_data_tables_cls(object) :
    file_name_start = 'chrom_'
    file_name_end = '_snp_data.h5'
    snp_data_table_name = 'snp_data'
    snp_bitpacked_allele_values_table_name = 'bitpacked_allele_values'
    def __init__(self, chrom) :
        self.chrom = chrom
        file_name = self.file_name_start + str(chrom) + self.file_name_end
        data_file_path = os.path.join(mod_dir, file_name)
        self.h5_file = tb.open_file(data_file_path, 'r')
        self.snp_data_table = getattr(self.h5_file.root, self.snp_data_table_name)
        self.snp_bitpacked_allele_values_table = getattr(self.h5_file.root, self.snp_bitpacked_allele_values_table_name)
                
    def close(self) :
        if self.h5_file is not None :
            self.h5_file.close()
            self.h5_file = None
            self.snp_data_table = None
            self.snp_bitpacked_allele_values_table = None
            
    def __del__(self) :
        self.close()
        

class snp_item_cls(object) :
    def __init__(self, item_data, item_allele_mask) :
        self.data = item_data
        self.allele_mask = item_allele_mask
        
class chrom_snp_item_factory_cls(object) :
    def __init__(self, chrom) :
        self.chrom = chrom
        self.data_tables = chrom_snp_data_tables_cls(chrom)

    def snp_item_objs(self, snp_indexes):
        snp_data = self.data_tables.snp_data_table[snp_indexes]
        bitpacked_alleles = self.data_tables.snp_bitpacked_allele_values_table[snp_indexes]
        allele_masks = np.unpackbits(bitpacked_alleles['bitpacked_values'], axis=1)
        snp_objs = []        
        for ind in xrange(snp_data.size) :
            snp_objs.append(snp_item_cls(snp_data[ind], allele_masks[ind]))
        return snp_objs

    def close(self) :
        if self.data_tables is not None :
            self.data_tables.close()
            self.data_tables = None
            
    def __del__(self) :
        self.close()
        
        
        
        
        
        
        
        
        
        
        