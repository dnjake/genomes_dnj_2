# -*- coding: utf-8 -*-
'''
I need to sort snps by pos before first and last pos selection
extending the last snp in a series can find earlier snps that meet the new criteria
'''

import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')
import numpy as np
from chrom_snp_series import snp_series_cls, single_snp_series_cls    
#from genome_data.sample_data_plus import sample_count_cls
#from genome_data.snp_data import snp_data_cls
from autosome_snp_data.chrom_snp_data_rdr import chrom_snp_data_tables_cls
#sc = sample_count_cls()
#sd_obj = snp_data_cls()

'''
may just need to modify this to take a chrom as input and create the right file and table names
'''
class chrom_snp_series_data_writer_cls(object) :
    total_allele_count = 5008
    data_dtype = np.dtype([('data_index', 'u4'), ('first_snp_index', 'u4'), ('chrom', 'u2'), ('first_pos', 'u4'),
                           ('last_pos', 'u4'), ('item_data_start', 'u4'), ('item_count', 'u2'), ('p90_allele_count', 'u2'), 
                           ('afr_obs_to_pred', 'f4'), ('afx_obs_to_pred', 'f4'), ('amr_obs_to_pred', 'f4'), 
                           ('eas_obs_to_pred', 'f4'), ('eur_obs_to_pred', 'f4'), ('sas_obs_to_pred', 'f4'), ('sax_obs_to_pred', 'f4')]  )
                           
    sample_values_dtype = np.dtype([('data_index', 'u2'), ('bitpacked_p90_allele_values', 'u1', (626,))])
    series_table_name_end = 'snp_series_data'
    series_items_table_name_end = 'series_items_snp_indexes'
    bitpacked_alleles_table_name_end = 'series_bitpacked_p90_allele_values'
    def __init__(self, chrom, series_rdr) :
        self.chrom = chrom
        self.table_writer = snp_series_data_table_writer_cls(chrom)
        self.snp_series = series_rdr.series_table[:]
        self.snp_series_items = series_rdr.series_items_table[:]        
        series_rdr.close()
        self.snp_data_rdr = chrom_snp_data_tables_cls(chrom)

    '''
    >>> ss_rdr.series_table.dtype
    dtype([('first_snp_index', '<u4'), ('item_data_start', '<u4'), ('item_count', '<u2')])
    >>> ss_rdr.series_item_table.dtype
    dtype([('snp_index', '<u4'), ('first_snp_index', '<u4')])
    >>> ss_rdr.single_snp_table.dtype
    dtype([('snp_index', '<u4')])    
    '''

    def write_out_data(self, data_index, item_data_start, item_count, series_obj) :
        series_obj.calc_p90_data()
        out_data = [data_index]
        out_data.extend(series_obj.get_series_range_data())
        out_data.extend([item_data_start, item_count])
        out_data.extend(series_obj.get_series_pop_data())
        self.out_data = out_data
        self.table_writer.append_series_data(tuple(out_data))
        self.table_writer.append_bitpacked_alleles(tuple([data_index, series_obj.get_bitpacked_p90_allele_mask()]))


    def write_series_data(self, data_index, index_snp_series) :
        first_snp_index, item_data_start, item_count = index_snp_series
        item_data_bound = item_data_start + item_count
        series_items = self.snp_series_items[item_data_start:item_data_bound]
        series_snp_indexes = series_items['snp_index']
        snp_data = self.snp_data_rdr.snp_data_table[series_snp_indexes]
        snp_bitpacked_alleles = self.snp_data_rdr.snp_bitpacked_allele_values_table[series_snp_indexes]
        chrom = self.chrom
        series_obj = snp_series_cls(chrom, snp_data, snp_bitpacked_alleles)
        self.write_out_data(data_index, item_data_start, item_count, series_obj)
        
    def write_single_snp_data(self, data_index, index_snp_series) :
        snp_index, item_data_start, item_count = index_snp_series
        snp_data = self.snp_data_rdr.snp_data_table[snp_index]
        snp_bitpacked_alleles = self.snp_data_rdr.snp_bitpacked_allele_values_table[snp_index]
        series_obj = single_snp_series_cls(self.chrom, snp_data, snp_bitpacked_alleles)
        self.write_out_data(data_index, item_data_start, item_count, series_obj)               

    def write_output(self) :
        for data_index, index_snp_series in enumerate(self.snp_series) :
            if index_snp_series['item_count'] == 1 :
                self.write_single_snp_data(data_index, index_snp_series)
            else :
                self.write_series_data(data_index, index_snp_series)
            if data_index % 10000 == 0 :
                print data_index            

    def do_write(self) :
        self.write_output()
        self.table_writer.close()
        self.snp_data_rdr.close()
        print 'done'


class snp_series_data_table_writer_cls(object) :

    file_name_end = '_snp_series_data.h5'
    series_data_table_name_end = '_snp_series_data'
    bitpacked_alleles_table_name_end = '_series_bitpacked_p90_allele_mask'
    
    data_dtype = np.dtype([('data_index', 'u4'), ('first_snp_index', 'u4'), ('chrom', 'u2'), ('first_pos', 'u4'),
                           ('last_pos', 'u4'), ('item_data_start', 'u4'), ('item_count', 'u2'), ('p90_allele_count', 'u2'), 
                           ('afr_obs_to_pred', 'f4'), ('afx_obs_to_pred', 'f4'), ('amr_obs_to_pred', 'f4'), 
                           ('eas_obs_to_pred', 'f4'), ('eur_obs_to_pred', 'f4'), ('sas_obs_to_pred', 'f4'), ('sax_obs_to_pred', 'f4')]  )                           
    
    allele_mask_dtype = np.dtype([('data_index', 'u4'), ('bitpacked_p90_allele_mask', 'u1', (626,))])


    def __init__(self, chrom) :
        self.chrom = chrom
        file_name = 'chrom_' + str(chrom) + self.file_name_end
        self.h5_out = tb.open_file(file_name, 'w', filters=filters)
        self.create_output_tables()

    def append_series_data(self, series_data) :
        self.series_data_table.append([series_data])
        
    def append_bitpacked_alleles(self, bitpacked_allele_mask) :
        self.bitpacked_alleles_table.append([bitpacked_allele_mask])
    
    def create_output_tables(self) :
        chrom = self.chrom
        chrom = str(chrom)
        h5_out = self.h5_out
        series_data_table_name = 'chrom_' + chrom + self.series_data_table_name_end
        self.series_data_table = h5_out.create_table('/', series_data_table_name, description=self.data_dtype)
        bitpacked_alleles_table_name = 'chrom_' + chrom + self.bitpacked_alleles_table_name_end
        self.bitpacked_alleles_table = h5_out.create_table('/', bitpacked_alleles_table_name, description=self.allele_mask_dtype)

    def close(self) :
        if self.h5_out is not None :
            self.h5_out.close()
            self.series_table = None
            self.bitpacked_alleles_table = None
            
    def __del__(self) :
        self.close()
        
        
        
        
        
        
        
        
        
        
        
        
        