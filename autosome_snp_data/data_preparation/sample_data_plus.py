# -*- coding: utf-8 -*-

'''
One addition for this module needs to be support to format the data in
a snp with sample measurements
'''



import os
import tables as tb
import numpy as np
from itertools import izip

s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)
data_file_path = os.path.join(mod_dir, 'phase3_sample_data.h5')

pop_names = {
    'CHB': 'China Bejing Han',
    'JPT': 'Japan Tokyo',
    'CHS': 'China Southern Han',
    'CDX': 'China Xishuangbanna Dai',
    'KHV': 'Vietnam Ho Chi Minh City Kinh',
    'CEU': 'European Utah',
    'TSI': 'Italy Toscani',
    'FIN': 'Finland',
    'GBR': 'Great Britain',
    'IBS': 'Spain Iberian',
    'YRI': 'Nigeria Ibadan Yoruba',
    'LWK': 'Kenya Webuye Luhya',
    'GWD': 'Gambia Western Divisions',
    'MSL': 'Sierra Leone Mende',
    'ESN': 'Nigeria Esan',
    'ASW': 'African Ancestry US South West',
    'ACB': 'African Caribbeans Barbados',
    'MXL': 'Mexican Ancestry US Los Angeles',
    'PUR': 'Puerto Rico',
    'CLM': 'Columbia Medellin',
    'PEL': 'Peru Lima',
    'GIH': 'Indian Gujarati US Houston',
    'PJL': 'Pakistan Lahore Punjabi',
    'BEB': 'Bangladesh Bengali',
    'STU': 'Sri Lankan Tamil UK',
    'ITU': 'Indian Telugu UK'        
}
    
super_pop_code_by_pop_code = {
    'CHB': 'EAS',
    'JPT': 'EAS',
    'CHS': 'EAS',
    'CDX': 'EAS',
    'KHV': 'EAS',
    'CEU': 'EUR',
    'TSI': 'EUR',
    'FIN': 'EUR',
    'GBR': 'EUR',
    'IBS': 'EUR',
    'YRI': 'AFR',
    'LWK': 'AFR',
    'GWD': 'AFR',
    'MSL': 'AFR',
    'ESN': 'AFR',
    'ASW': 'AFX',
    'ACB': 'AFX',
    'MXL': 'AMR',
    'PUR': 'AMR',
    'CLM': 'AMR',
    'PEL': 'AMR',
    'GIH': 'SAX',
    'PJL': 'SAS',
    'BEB': 'SAS',
    'STU': 'SAX',
    'ITU': 'SAX'        
}

super_pop_plus_names = {
    'AFR': 'Africa',
    'AFX': 'African Ancestry',
    'AMR': 'America',
    'EAS': 'East Asian',
    'EUR': 'European',
    'SAS': 'South Asian',
    'SAX': 'South Asian Ancestry'
}

super_pop_plus_nums_in = [
    ('AFR', 0), 
    ('AFX', 1),
    ('AMR', 2), 
    ('EAS', 3), 
    ('EUR', 4), 
    ('SAS', 5),
    ('SAX', 6)
]

super_pop_plus_dt = np.dtype([('super_pop', 'S3'),('super_pop_plus_num', 'u2')])

pop_num_map_in = [
       ('ACB', 0, 1), ('ASW', 1, 1), ('BEB', 2, 5), ('CDX', 3, 3), ('CEU', 4, 4),
       ('CHB', 5, 3), ('CHS', 6, 3), ('CLM', 7, 2), ('ESN', 8, 0), ('FIN', 9, 4),
       ('GBR', 10, 4), ('GIH', 11, 6), ('GWD', 12, 0), ('IBS', 13, 4), ('ITU', 14, 6),
       ('JPT', 15, 3), ('KHV', 16, 3), ('LWK', 17, 0), ('MSL', 18, 0), ('MXL', 19, 2),
       ('PEL', 20, 2), ('PJL', 21, 5), ('PUR', 22, 2), ('STU', 23, 6), ('TSI', 24, 4),
       ('YRI', 25, 0)
   ]
   
pop_map_dt = np.dtype([('pop', 'S3'), ('pop_num', 'u2'), ('super_pop_plus_num', 'u2')])


def super_pop_plus_sample_counts_def(pop_num_map, pop_stats) :
    num_super_pops = pop_num_map['super_pop_plus_num'].max() + 1
    samples_in_super_pop = np.zeros(num_super_pops, 'f4')
    super_pop_plus_nums = pop_num_map['super_pop_plus_num']
    for pop_num, pop, samples_in_pop, samples_freq in pop_stats :
        super_pop_num = super_pop_plus_nums[pop_num]
        samples_in_super_pop[super_pop_num] += samples_in_pop
    return samples_in_super_pop

class sample_count_cls(object) :
    pop_anal_result_descr = [('pop', 'S3'), ('allele_count', '<f4'), ('freq', '<f4'), ('obs_to_pred', '<f4')]
    pop_anal_result_dtype = np.dtype(pop_anal_result_descr)
    h_in = tb.open_file(data_file_path, 'r')
    root = h_in.root
    # Provides the thousand genome phase 3 identification for each sample in the data
    # A sample is a person.  This table associates each sample with and index from 0 to 2503
    samples = root.integrated_call_samples_panel[:]
    # The pop codes give the 3 letter code for the thousand genome population
    # The index of the code in this table is used to identify the population in sample tables
    pop_codes = root.pop_codes[:]
    # Equivalent data for super pops
    super_pop_codes = root.super_pop_codes[:]
    # pop stats gives count of number of samples in a population
    # A sample is a person in this table
    # There are 2504 samples in the thousand genome phase 3 data
    pop_stats = root.pop_stats[:]
    alleles_in_pop = 2.0*(pop_stats['samples_in_pop'].astype('f4'))
    # super pop stats provides equivalent counts for regions
    super_pop_stats = root.super_pop_stats[:]
    # Table of phases by offset in phase 3 thousand genome data
    # Phases are associated with an index from 0 to 5007
    # Index of a phase is 2 times the index of the sample 
    # plus 0 or 1 for the two sample phases
    # Table gives pop and super pop indexes for the phase
    phased_allele_calls = root.phased_sample_calls[:]
    h_in.close()
    super_pop_plus_nums = np.array(super_pop_plus_nums_in, super_pop_plus_dt)
    pop_num_map = np.array(pop_num_map_in, pop_map_dt)
    alleles_in_super_pop_plus = 2.0*super_pop_plus_sample_counts_def(pop_num_map, pop_stats)
    total_allele_count = alleles_in_super_pop_plus.sum()
    super_pop_plus_allele_freqs = alleles_in_super_pop_plus/total_allele_count
    pop_count_descr = []
    pop_names = pop_codes['code']
    for pop in pop_names :
        pop = pop.lower()
        pop_count_descr.append((pop + '_count', 'u2'))
        
    pop_count_dtype = np.dtype(pop_count_descr)
    pop_count_names = pop_count_dtype.names
    super_pop_plus_obs_to_freq_descr = []
    super_pop_plus_names = super_pop_plus_nums['super_pop']
    for super_pop in super_pop_plus_names :
        super_pop = super_pop.lower()
        super_pop_plus_obs_to_freq_descr.append((super_pop + '_obs_to_pred', 'f4'))
    super_pop_plus_obs_to_freq_dtype = np.dtype(super_pop_plus_obs_to_freq_descr)
    super_pop_plus_obs_to_freq_names = super_pop_plus_obs_to_freq_dtype.names
    item_data_dtype = np.dtype([('name', 'S3'), ('allele_count', 'u2'), ('obs_to_pred', 'f4') ])
    pop_data_size = len(pop_count_names)
    #total_allele_count = 5008.0
    pop_allele_freqs = alleles_in_pop/total_allele_count
    #super_pop_plus_allele_freqs = alleles_in_super_pop_plus/total_allele_count
    def __init__(self) :
        pass
    
    def pop_allele_counts_from_allele_indexes(self, allele_indexes) :
        num_pops = self.pop_codes.size
        pop_allele_counts = np.zeros(num_pops, 'f4')
        allele_index_pop_nums = self.phased_allele_calls['pop_num']
        for ind in allele_indexes :
            pop_allele_counts[allele_index_pop_nums[ind]] += 1
        return pop_allele_counts
        
    def super_pop_plus_from_pop_allele_counts(self, pop_allele_counts) :
        num_super_pop_plus = self.super_pop_plus_nums.size
        counts = np.zeros(num_super_pop_plus, 'f4')
        pop_super_pop_plus_nums = self.pop_num_map['super_pop_plus_num']
        for ind, count in enumerate(pop_allele_counts) :
            super_pop_plus_num = pop_super_pop_plus_nums[ind]
            counts[super_pop_plus_num] += count
        return counts

    def super_pop_plus_and_pop_stats(self, allele_indexes) :
        total_allele_count = float(allele_indexes.size)
        pop_counts = self.pop_allele_counts_from_allele_indexes(allele_indexes)
        # for super_pop_plus_counts I want the freq_to_pred
        super_pop_plus_counts = self.super_pop_plus_from_pop_allele_counts(pop_counts)
        super_pop_plus_pred_counts = total_allele_count*self.super_pop_plus_allele_freqs
        super_pop_plus_obs_to_pred = super_pop_plus_counts/super_pop_plus_pred_counts
        return super_pop_plus_obs_to_pred, pop_counts

    def pop_data(self, snp_allele_data_item) :
        pop_count_data = []
        for pop_name, pop_count_name in izip(self.pop_names, self.pop_count_names) :
            pop_count = snp_allele_data_item[pop_count_name]
            pop_count_data.append((pop_name, pop_count, 0.0))
        pop_count_data = np.array(pop_count_data, self.item_data_dtype)
        snp_alleles = float(pop_count_data['allele_count'].sum())
        snp_freq = snp_alleles/self.total_allele_count
        pop_pred_alleles = snp_freq*self.alleles_in_pop
        pop_count_data['obs_to_pred'] = (pop_count_data['allele_count'].astype('f4'))/pop_pred_alleles
        return pop_count_data


    def super_pop_plus_data(self, snp_allele_data_item) :
        obs_to_pred_data = []
        for super_pop_plus_name , obs_to_pred_name in izip(self.super_pop_plus_names, self.super_pop_plus_obs_to_freq_names ) :
            obs_to_pred = snp_allele_data_item[obs_to_pred_name]
            obs_to_pred_data.append((super_pop_plus_name, 0, obs_to_pred))
        obs_to_pred_data = np.array(obs_to_pred_data, self.item_data_dtype)
        snp_alleles = snp_allele_data_item['all_count']
        snp_freq = snp_alleles/self.total_allele_count
        pred_alleles = snp_freq*self.alleles_in_super_pop_plus
        obs_alleles = obs_to_pred_data['obs_to_pred']*pred_alleles
        obs_alleles += 0.5
        obs_to_pred_data['allele_count'] = obs_alleles
        return obs_to_pred_data
        
    def super_pop_plus_and_pop_data(self, snp_allele_data_item) :
        pd = self.pop_data(snp_allele_data_item)
        sppd = self.super_pop_plus_data(snp_allele_data_item)
        return sppd, pd

    def analyze_countries(self, allele_indexes) :
        num_pops = self.pop_codes.size
        #allele_index_pop_nums = self.phased_allele_calls['pop_num']
        anal_a = np.zeros(num_pops, self.pop_anal_result_dtype)
        anal_a['pop'] = self.pop_codes['code']
        anal_a['allele_count'] = self.pop_allele_counts_from_allele_indexes(allele_indexes)
        anal_a['freq'] = anal_a['allele_count'] / self.alleles_in_pop
        exp_alleles_in_pop = float(allele_indexes.size)*self.pop_allele_freqs
        anal_a['obs_to_pred'] = anal_a['allele_count'] / exp_alleles_in_pop
        return anal_a

    def analyze_regions_from_pops(self, allele_indexes, pop_allele_counts) :
        num_super_pop_plus = self.super_pop_plus_nums.size
        #sample_index_super_pop_nums = self.phased_sample_calls['super_pop_num']
        anal_a = np.zeros(num_super_pop_plus, np.dtype(self.pop_anal_result_dtype))
        anal_a['pop'] = self.super_pop_plus_names
        # The allele count is the number of alleles in each of the regional populations
        anal_a['allele_count'] = self.super_pop_plus_from_pop_allele_counts(pop_allele_counts)
        # The freq is just the allele count dividied by the number of alleles in the population
        anal_a['freq'] = anal_a['allele_count'] / self.alleles_in_super_pop_plus
        # calculates the expected population alleles from the total number of alleles and the pop freq
        exp_alleles_in_super_pop_plus = float(allele_indexes.size)*self.super_pop_plus_allele_freqs
        anal_a['obs_to_pred'] = anal_a['allele_count'] / exp_alleles_in_super_pop_plus
        return anal_a

    def analyze_regions(self, allele_indexes) :
        pop_allele_counts = self.pop_allele_counts_from_allele_indexes(allele_indexes)
        return self.analyze_regions_from_pops(allele_indexes, pop_allele_counts)
        
    def analyze_regions_and_countries(self, allele_indexes) :
        pop_data = self.analyze_countries(allele_indexes)
        super_pop_plus_data = self.analyze_regions_from_pops(allele_indexes, pop_data['allele_count'])
        return super_pop_plus_data, pop_data






