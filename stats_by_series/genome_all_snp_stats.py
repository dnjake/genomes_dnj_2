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


stats_file_name = 'genome_all_snp_stats.h5'
stats_table_name = 'genome_all_snp_stats'
stats_file_path = os.path.join(mod_dir, stats_file_name)

stats_dtype = np.dtype([('chrom', 'u2'), ('snp_1', 'u4'), ('snp_2', 'u4'), ('snp_3', 'u4'), ('non_series_snps', 'u4'),
                        ('series_count', 'u4'), ('series_snps', 'u4'), ('mean_snps_per_series', 'f4'),
                        ('median_snps_per_series', 'u2'), ('max_snps_per_series', 'u2'),
                        ('mean_samples_per_series', 'f4'), ('median_samples_per_series', 'u2'),
                        ('mean_series_length', 'f4'), ('median_series_length', 'u4'), ('max_series_length', 'u4'),
                        ('no_afr_series_count', 'u4'), ('low_afr_series_count', 'u4'), ('some_afr_series_count', 'u4'),
                        ('high_afr_series_count', 'u4'), ('all_afr_series_count', 'u4')])


class chrom_snp_series_stats_cls(object) :
    def __init__(self, chrom) :
        self.chrom = chrom        
        rdr = data_mod.chrom_snp_series_data_tables_cls(chrom)
        self.series_data = rdr.series_data_table[:]
        rdr.close()

    def calc_series_stats(self, stats) :    
        m = self.series_data['item_count'] >= 4
        in_series_data = self.series_data[m]
        stats.append(in_series_data.size)
        stats.extend(ssh.snps_count_mean_median_max(in_series_data))
        stats.extend(ssh.samples_mean_median(in_series_data))
        stats.extend(ssh.lengths_mean_median_max(in_series_data))
        stats.extend(ssh.series_counts_by_pop_type(in_series_data))
    
    def calc_non_series_snps(self, stats) :                     
        non_series_snps = 0
        for count in range(1,4) :
            m = self.series_data['item_count'] == count
            count_series = m.sum()
            stats.append(count_series)
            non_series_snps += count_series*count
        stats.append(non_series_snps)
        
    def calc_stats(self, stats) :
        self.calc_non_series_snps(stats)
        self.calc_series_stats(stats)
        
    def do_stats(self, stats_table) :
        stats = [self.chrom]
        self.calc_stats(stats)
        stats_table.append([tuple(stats)])

class genome_snp_series_stats_cls(object) :
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
            chrom_obj = chrom_snp_series_stats_cls(chrom)
            chrom_obj.do_stats(self.stats_table)
        self.close_h5()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
