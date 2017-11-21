# -*- coding: utf-8 -*-

import numpy as np
import tables as tb
import os

from ..autosome_snp_data import chrom_snp_series_data_rdr as data_mod
from ..stats_by_series import series_stats_helpers as ssh

s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

def h5_filters() :
    return tb.Filters(complib='zlib', complevel=5)


stats_file_name = 'genome_series_event_stats.h5'
stats_table_name = 'genome_series_event_stats'
stats_file_path = os.path.join(mod_dir, stats_file_name)

stats_dtype = np.dtype([('chrom', 'u2'), ('event_first_pos', 'u4'), ('event_last_pos', 'u4'), 
                        ('series_min_first_pos', 'u4'), ('series_max_last_pos', 'u4'), ('mean_series_length', 'f4'),
                        ('samples_in_series', 'u2'), 
                        ('series_count', 'u2'), ('sample_weighted_length', 'u8'),
                        ('snps_in_series', 'u4'), ('mean_series_snps', 'f4'), ('sample_weighted_snps', 'u8'),
                        ('no_afr_series_count', 'u2'), ('no_afr_sample_weighted_length', 'u8'),
                        ('low_afr_series_count', 'u2'), ('low_afr_sample_weighted_length', 'u8'),
                        ('some_afr_series_count', 'u2'), ('some_afr_sample_weighted_length', 'u8'),
                        ('high_afr_series_count', 'u2'), ('high_afr_sample_weighted_length', 'u8'), 
                        ('all_afr_series_count', 'u2'), ('all_afr_sample_weighted_length', 'u8'),
                        ])
                        
                        
base_non_pos_data_null = [0, 0, 0.0,
                          0,
                          0, 0,  
                          0, 0.0, 0, 
                          0, 0, 
                          0, 0, 
                          0, 0, 
                          0, 0, 
                          0, 0] 

class chrom_series_event_stats_cls(object) :    
    non_pos_stats_null = []
    non_pos_stats_null.extend(base_non_pos_data_null)
    def __init__(self, chrom) :
        self.chrom = chrom        
        rdr = data_mod.chrom_snp_series_data_tables_cls(chrom)
        series_data = rdr.series_data_table[:]
        m = series_data['item_count'] >= 4
        self.series_data = series_data[m]
        bp_alleles = rdr.bitpacked_alleles_table[:]
        bp_alleles = bp_alleles[m]
        rdr.close()
        self.series_p90_allele_masks = np.unpackbits(bp_alleles['bitpacked_p90_allele_mask'], axis=1)
        self.first_pos = self.series_data['first_pos']
        self.last_pos = self.series_data['last_pos']
        self.active_series = {}
        self.last_pos_sa = self.series_data['last_pos'].argsort()
        self.next_first_pos_array_index = 0
        self.next_last_pos_sa_index = 0

    def do_pop_types(self, active_series_data, event_stats) :
        pop_type_series_data = ssh.series_data_no_low_some_high_all_afr(active_series_data)
        for sd in pop_type_series_data :
            if sd.size == 0 :
                event_stats.extend([0, 0])
            else :
                event_stats.extend(ssh.series_count_sample_weighted_length(sd))
        

    def process_event(self, event_stats) :
        active_series_array_indexes = self.active_series.keys()
        assert len(active_series_array_indexes) == self.active_series_count
        active_series_data = self.series_data[active_series_array_indexes]
        active_allele_masks = self.series_p90_allele_masks[active_series_array_indexes]
        event_stats.extend(ssh.pos_min_max_length_mean(active_series_data))
        event_stats.append(ssh.samples_in_series(active_allele_masks))
        event_stats.extend(ssh.series_count_sample_weighted_length(active_series_data))
        event_stats.extend(ssh.snps_count_mean_sample_weighted(active_series_data))
        self.do_pop_types(active_series_data, event_stats)

    def zero_stats(self, first_pos, last_pos) :
        stats_v = [self.chrom, first_pos, last_pos]
        stats_v.extend(self.non_pos_stats_null)
        stats_v = tuple(stats_v)
        #print stats_v
        self.stats_table.append([stats_v])

    def is_first_pos_next(self) :
        next_last_pos = self.last_pos[self.last_pos_sa[self.next_last_pos_sa_index]]
        if self.next_first_pos_array_index == self.first_pos.size :
            return False, next_last_pos
        next_first_pos = self.first_pos[self.next_first_pos_array_index]
        if next_first_pos < next_last_pos :
            return True, next_first_pos
        else :
            return False, next_last_pos

    def set_initial_next(self) :
        self.active_series_count = 0
        event_first_pos = 0
        self.event_last_pos = self.first_pos[0]
        self.next_first_pos_array_index = 0
        self.next_last_pos_sa_index = 0
        self.next_event_is_series_start = True
        self.zero_stats(event_first_pos, self.event_last_pos)


    def find_next_pos(self) :
        event_first_pos = self.event_last_pos
        event_is_series_start = self.next_event_is_series_start
        #print event_first_pos, event_is_series_start
        if event_is_series_start :
            event_series_array_index = self.next_first_pos_array_index
            #print event_series_array_index, event_first_pos
            assert event_first_pos == self.first_pos[event_series_array_index]
            self.next_first_pos_array_index += 1
        else :
            event_series_array_index = self.last_pos_sa[self.next_last_pos_sa_index]
            assert event_first_pos == self.last_pos[event_series_array_index]
            self.next_last_pos_sa_index += 1
        #print event_is_series_start, event_series_array_index            
        if event_is_series_start :
            #print event_series_array_index
            self.active_series[event_series_array_index] = True
            self.active_series_count += 1
        else :
            #print 'deleted', event_series_array_index
            del self.active_series[event_series_array_index]
            self.active_series_count -= 1
        self.next_event_is_series_start, self.event_last_pos = self.is_first_pos_next()
        return event_first_pos, self.event_last_pos

    def process_chrom(self, stats_table) :
        self.stats_table = stats_table
        self.set_initial_next()
        count = 0
        while self.next_last_pos_sa_index < (self.last_pos_sa.size-1) :
            event_first_pos, event_last_pos = self.find_next_pos()
            if self.active_series_count == 0 :
                self.zero_stats(event_first_pos, event_last_pos)
            else :
                event_stats = [self.chrom, event_first_pos, event_last_pos]
                self.process_event(event_stats)
                #print event_stats
                self.stats_table.append([tuple(event_stats)])
            if count % 1000 == 0 :
                print count
            count += 1



class genome_series_event_stats_writer_cls(object) :    
    def open_h5(self) :
        self.h5 = tb.open_file(stats_file_path, 'w', filters=h5_filters())
        self.stats_table = self.h5.create_table('/', stats_table_name, description=stats_dtype)
        
    def close_h5(self) :
        self.h5.close()
        self.h5 = None
        self.stats_table = None

    def do_all_chrom_stats(self) :
        self.open_h5()
        for chrom in range(1, 23) :
            print 'processing', chrom
            self.chrom_handler = chrom_series_event_stats_cls(chrom)
            self.chrom_handler.process_chrom(self.stats_table)            
        self.close_h5()
















































                        