# -*- coding: utf-8 -*-


import numpy as np
import tables as tb
import os

from ..autosome_snp_data import chrom_snp_series_data_rdr as data_mod
from . import series_stats_helpers as ssh

s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

def h5_filters() :
    return tb.Filters(complib='zlib', complevel=5)


stats_file_name = 'genome_by_series_stats.h5'
stats_table_name = 'genome_by_series_stats'
stats_file_path = os.path.join(mod_dir, stats_file_name)

stats_dtype = np.dtype([('chrom', 'u2'), ('length', 'u4'), ('snp_count', 'u2'),
                        ('sample_count', 'u2'), ('pop_type', 'u2')])
pop_types = {'no_afr': 1, 'low_afr': 2, 'some_afr': 3, 'high_afr':4, 'all_afr': 5}

class chrom_by_series_stats_cls(object) :
    min_snps_per_series = 4
    def __init__(self, chrom) :
        self.chrom = chrom        
        rdr = data_mod.chrom_snp_series_data_tables_cls(chrom)
        series_data = rdr.series_data_table[:]
        rdr.close()
        m = series_data['item_count'] >= self.min_snps_per_series
        self.series_data = series_data[m]

    def find_pop_types(self, data) :
        dpt = data['pop_type']
        m_no_afr = ssh.mask_no_afr(self.series_data)
        dpt[m_no_afr] = pop_types['no_afr']
        m_low_afr = ssh.mask_low_afr(self.series_data)
        dpt[m_low_afr] = pop_types['low_afr']
        m_some_afr = ssh.mask_some_afr(self.series_data)
        dpt[m_some_afr] = pop_types['some_afr']
        m_all_afr = ssh.mask_all_afr(self.series_data)
        dpt[m_all_afr] = pop_types['all_afr']
        m_high_afr = ssh.mask_high_afr(self.series_data, m_all_afr)
        dpt[m_high_afr] = pop_types['high_afr']

    def calc_stats(self) :
        data = np.zeros(self.series_data.size, stats_dtype)
        data['chrom'] = self.chrom
        data['length'] = self.series_data['last_pos'] - self.series_data['first_pos'] + 1
        data['snp_count'] = self.series_data['item_count']
        data['sample_count'] = self.series_data['p90_allele_count']
        self.find_pop_types(data)
        return data

    def do_stats(self, stats_table) :
        stats_table.append(self.calc_stats())

class genome_by_series_stats_cls(object) :
    def open_h5(self) :
        self.h5 = tb.open_file(stats_file_path, 'w', filters=h5_filters())
        self.stats_table = self.h5.create_table('/', stats_table_name, description=stats_dtype)
        
    def close_h5(self) :
        self.h5.close()
        self.h5 = None
        self.stats_table = None

    def do_stats(self) :
        self.open_h5()
        for chrom in range(1,23) :
            chrom_obj = chrom_by_series_stats_cls(chrom)
            chrom_obj.do_stats(self.stats_table)
        self.close_h5()
