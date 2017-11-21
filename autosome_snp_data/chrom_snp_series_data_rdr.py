# -*- coding: utf-8 -*-

import tables as tb
import numpy as np
import os


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)


class chrom_snp_series_data_tables_cls(object) :
    series_data_dtype = np.dtype([('data_index', '<u4'), ('first_snp_index', '<u4'), ('chrom', '<u2'), ('first_pos', '<u4'), 
                                  ('last_pos', '<u4'), ('item_data_start', '<u4'), ('item_count', '<u2'),
                                  ('p90_allele_count', '<u2'), ('afr_obs_to_pred', '<f4'), ('afx_obs_to_pred', '<f4'), 
                                  ('amr_obs_to_pred', '<f4'), ('eas_obs_to_pred', '<f4'), ('eur_obs_to_pred', '<f4'), 
                                  ('sas_obs_to_pred', '<f4'), ('sax_obs_to_pred', '<f4')])
                                  
    bitpacked_alleles_dtype = np.dtype([('data_index', '<u2'), ('bitpacked_p90_allele_mask', 'u1', (626,))])  
    file_name_end = '_snp_series_data.h5'
    series_data_table_name_end = '_snp_series_data'
    bitpacked_alleles_table_name_end = '_series_bitpacked_p90_allele_mask'
    def __init__(self, chrom) :
        self.chrom = chrom
        table_name_start = 'chrom_' + str(chrom)
        series_data_table_name = table_name_start + self.series_data_table_name_end
        bitpacked_alleles_table_name = table_name_start + self.bitpacked_alleles_table_name_end
        file_name = 'chrom_' + str(chrom) + self.file_name_end
        file_path = os.path.join(mod_dir, file_name)
        self.h5 = tb.open_file(file_path, 'r')
        self.series_data_table = getattr(self.h5.root, series_data_table_name)
        self.bitpacked_alleles_table = getattr(self.h5.root, bitpacked_alleles_table_name)
                
    def from_data_indexes(self, data_indexes) :
        series_data = self.series_data_table[data_indexes]
        bp_alleles = self.bitpacked_alleles_table[data_indexes]
        allele_masks = np.unpackbits(bp_alleles['bitpacked_p90_allele_mask'], axis=1)
        return series_data, allele_masks

    def close(self) :
        if self.h5 is not None :
            self.h5.close()
            self.h5 = None
            self.series_data_table = None
            self.single_snp_table = None
            
    def __del__(self) :
        self.close()
        

        
class chrom_pos_range_series_rdr_cls(object) :
    def __init__(self, chrom, range_start_pos, range_bound_pos, min_allele_count=None) :
        self.chrom = chrom
        self.range_start_pos = range_start_pos
        self.range_bound_pos = range_bound_pos
        self.min_allele_count = min_allele_count
        gcond = '(last_pos > ' + str(range_start_pos) + ')'
        lcond = '(first_pos < ' + str(range_bound_pos) + ')'
        self.range_condition = gcond + ' & ' + lcond
        
    def read_range(self) :
        tables = chrom_snp_series_data_tables_cls(self.chrom)        
        self.series_data = tables.series_data_table.read_where(self.range_condition)
        if self.min_allele_count is not None :
            m = self.series_data['p90_allele_count'] >= self.min_allele_count
            self.series_data = self.series_data[m]
        data_indexes = self.series_data['data_index']
        bp_alleles = tables.bitpacked_alleles_table[data_indexes]        
        self.allele_masks = np.unpackbits(bp_alleles['bitpacked_p90_allele_mask'], axis=1)
        tables.close()


class chrom_series_in_pos_interval_rdr_cls(object) :
    def __init__(self, chrom, interval_start_pos, interval_bound_pos, min_allele_count=None) :
        self.chrom = chrom
        self.interval_start_pos = interval_start_pos
        self.interval_bound_pos = interval_bound_pos
        self.min_allele_count = min_allele_count
        gcond = '(first_pos >= ' + str(interval_start_pos) + ')'
        lcond = '(last_pos <= ' + str(interval_bound_pos) + ')'
        self.range_condition = gcond + ' & ' + lcond
        
    def read_interval_series(self) :
        tables = chrom_snp_series_data_tables_cls(self.chrom)        
        self.series_data = tables.series_data_table.read_where(self.range_condition)
        if self.min_allele_count is not None :
            m = self.series_data['p90_allele_count'] >= self.min_allele_count
            self.series_data = self.series_data[m]
        data_indexes = self.series_data['data_index']
        bp_alleles = tables.bitpacked_alleles_table[data_indexes]        
        self.allele_masks = np.unpackbits(bp_alleles['bitpacked_p90_allele_mask'], axis=1)
        tables.close()
 
























           