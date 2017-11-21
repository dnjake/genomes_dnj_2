# -*- coding: utf-8 -*-

import numpy as np
#from genome_data.sample_data_plus import sample_count_cls
from genomes_dnj.autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
'''
output data for snp series  The series object will create everything but the data index
maybe I will pass in the data index so the object can create the whole tuple
    data_dtype = np.dtype([('data_index', 'u4'), ('first_snp_index', 'u4'), ('chrom', 'u2'), ('first_pos', 'u4'),
                           ('last_pos', 'u4'), ('item_data_start', 'u4'), ('item_count', 'u2'), ('p90_allele_count', 'u2'), 
                           ('afr_obs_to_pred', 'f4'), ('afx_obs_to_pred', 'f4'), ('amr_obs_to_pred', 'f4'), 
                           ('eas_obs_to_pred', 'f4'), ('eur_obs_to_pred', 'f4'), ('sas_obs_to_pred', 'f4'), ('sax_obs_to_pred', 'f4')]  )
'''

class snp_series_cls(object) :
    #pda = sample_count_cls()
    pda = country_region_alleles_cls()
    total_allele_count = 5008
    struct_data_descr = [ ('chrom', 'u2'), ('first_snp_index', 'i4'), ('series_snp_count', 'u2'),  ('p90_allele_count', 'u2'),  
                         ('afr_obs_to_pred', 'f4'), ('afx_obs_to_pred', 'f4'), ('amr_obs_to_pred', 'f4'), 
                         ('eas_obs_to_pred', 'f4'), ('eur_obs_to_pred', 'f4'), 
                         ('sas_obs_to_pred', 'f4'), ('sax_obs_to_pred', 'f4'), ('snp_series_obj', 'O') 
    ]
    def __init__(self, chrom, series_snp_data, series_bitpacked_alleles) :
        self.chrom = chrom
        self.series_snp_data = series_snp_data
        self.series_bitpacked_alleles = series_bitpacked_alleles
        self.series_snp_indexes = self.series_snp_data['snp_index']
        self.snp_series_pos = self.series_snp_data['pos']
        
    #('first_snp_index', 'u4'), ('chrom', 'u2'), ('first_pos', 'u4'), ('last_pos', 'u4')
    def get_series_range_data(self) :
        first_snp_index = self.series_snp_indexes[0]
        first_pos = self.snp_series_pos[0]
        last_pos  = self.snp_series_pos[-1]
        return [first_snp_index, self.chrom, first_pos, last_pos]
        
        
    def calc_p90_data(self) :
        series_alleles_that_expr_snp = np.unpackbits(self.series_bitpacked_alleles['bitpacked_values'], axis=1)
        series_snps_per_allele = np.sum(series_alleles_that_expr_snp, axis=0)
        snp_count = self.series_snp_indexes.size
        p90_snp_count = 0.9*float(snp_count) + 0.5
        p90_snp_count = int(p90_snp_count)
        p90_mask = series_snps_per_allele >= p90_snp_count
        self.p90_allele_indexes = np.where(p90_mask)[0]
        self.p90_super_pop_plus_data = self.pda.analyze_regions(self.p90_allele_indexes)
        

    def get_series_pop_data(self) :
        out_data = [self.p90_allele_indexes.size]
        out_data.extend(self.p90_super_pop_plus_data['obs_to_pred'].tolist())
        return out_data

    def get_bitpacked_p90_allele_mask(self) :
        p90_allele_mask = np.zeros(self.total_allele_count, 'u1')
        p90_allele_mask[self.p90_allele_indexes] = 1
        return np.packbits(p90_allele_mask)
        

class single_snp_series_cls(snp_series_cls) :
    def __init__(self, chrom, snp_data, snp_bitpacked_alleles) :
        self.chrom = chrom
        self.snp_data = snp_data
        self.snp_index = self.snp_data['snp_index']
        self.snp_bitpacked_alleles = snp_bitpacked_alleles
        
    def get_series_range_data(self) :
        first_snp_index = self.snp_index
        first_pos = self.snp_data['pos']
        last_pos = first_pos
        return [first_snp_index, self.chrom, first_pos, last_pos]

    def calc_p90_data(self) :
        alleles_that_expr_snp = np.unpackbits(self.snp_bitpacked_alleles['bitpacked_values'])
        allele_indexes = np.where(alleles_that_expr_snp)[0]
        self.p90_allele_indexes = allele_indexes
        self.p90_super_pop_plus_data = self.pda.analyze_regions(self.p90_allele_indexes)





