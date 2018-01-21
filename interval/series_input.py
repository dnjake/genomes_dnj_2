# -*- coding: utf-8 -*-

import numpy as np
from ..stats_by_pos import zero_series_pos as zps
from ..autosome_snp_data.chrom_snp_series_data_rdr import chrom_series_in_pos_interval_rdr_cls
from . import series_anal as sa

class series_input_selection_cls(object) :
    def __init__(self, chrom, series_data, sample_data, first_pos=None, last_pos=None) :
        self.chrom = chrom
        self.series_data = series_data
        self.sample_data = sample_data
        self.allele_masks = self.sample_data['allele_mask']
        self.in_first_pos = first_pos
        self.in_last_pos = last_pos
        self.series_first_pos = self.series_data['first_pos']
        self.series_last_pos = self.series_data['last_pos']
        
    def full_input_interval(self) :
        series_sample_data = {'series_data': self.series_data, 'sample_data': self.sample_data}
        return sa.series_data_anal_cls(self.chrom, series_sample_data)
        
    def overlap_pos(self, pos) :
        m = self.series_first_pos <= pos
        m = np.logical_and(m, self.series_last_pos >= pos)
        series_sample_data = {'series_data': self.series_data[m], 'sample_data': self.sample_data[m]}
        return sa.series_data_anal_cls(self.chrom, series_sample_data)
        
    def overlap_interval(self, interval_range) :
        interval_first_pos, interval_last_pos = interval_range
        m = self.series_first_pos <= interval_last_pos
        m = np.logical_and(m, self.series_last_pos >= interval_first_pos)
        series_sample_data = {'series_data': self.series_data[m], 'sample_data': self.sample_data[m]}
        out_obj = sa.series_data_anal_cls(self.chrom, series_sample_data, 
                                          selection_range=(interval_first_pos, interval_last_pos))
        return out_obj
                
    def within_interval(self, interval_first_pos, interval_last_pos) :
        m = self.series_first_pos >= interval_first_pos
        m = np.logical_and(m, self.series_first_pos < interval_last_pos)
        m = np.logical_and(m, self.series_last_pos > interval_first_pos)
        m = np.logical_and(m, self.series_last_pos <= interval_last_pos)
        series_sample_data = {'series_data': self.series_data[m], 'sample_data': self.sample_data[m]}
        return sa.series_data_anal_cls(self.chrom, series_sample_data)

    def overlap_in_first_pos(self) :
        return self.overlap_pos(self.in_first_pos)
        
    def overlap_in_first_last_pos(self) :
        return self.overlap_interval((self.in_first_pos, self.in_last_pos))
        
    def within_in_first_last_pos(self) :
        return self.within_interval(self.in_first_pos, self.in_last_pos)
        
class series_data_rdr_cls(object) :
    sample_data_dtype = sa.series_data_anal_cls.sample_data_dtype
    min_series_snp_count = 4
    def __init__(self, chrom) :
        self.chrom = chrom
        zps_rdr = zps.chrom_zero_pos_cls(chrom)
        self.chrom_zero_pos = zps_rdr.chrom_zero_pos
        self.zero_pos = self.chrom_zero_pos['pos']

    def find_interval_first_pos(self, first_pos) :
        ind_zero_pos = self.zero_pos.searchsorted(first_pos)
        interval_first_pos = self.zero_pos[ind_zero_pos]
        if first_pos < interval_first_pos :
            ind_zero_pos -= 1
            interval_first_pos = self.zero_pos[ind_zero_pos]
        return interval_first_pos
        
    def find_interval_last_pos(self, last_pos) :
        ind_zero_pos = self.zero_pos.searchsorted(last_pos)
        interval_last_pos = self.zero_pos[ind_zero_pos]
        return interval_last_pos

    def read_zero_interval_data(self, interval_first_pos, interval_last_pos) :
        '''
        A series that ends at the last pos will be included in the read
        '''
        rdr = chrom_series_in_pos_interval_rdr_cls(self.chrom, interval_first_pos, interval_last_pos)
        rdr.read_interval_series()
        m = rdr.series_data['item_count'] >= self.min_series_snp_count
        series_data = rdr.series_data[m]
        sample_data = np.zeros(series_data.size, self.sample_data_dtype)
        sample_data['data_index'] = series_data['data_index']
        sample_data['allele_mask'] = rdr.allele_masks[m]
        return series_data, sample_data
    
    def read_series_data_for_interval(self, first_pos, last_pos=None) :
        interval_first_pos = self.find_interval_first_pos(first_pos)
        if last_pos is None :
            last_pos = first_pos
        interval_last_pos = self.find_interval_last_pos(last_pos)
        return self.read_zero_interval_data(interval_first_pos, interval_last_pos)
    
    def input_selection_obj(self, first_pos, last_pos=None) :
        series_data, sample_data = self.read_series_data_for_interval(first_pos, last_pos)
        return series_input_selection_cls(self.chrom, series_data, sample_data, first_pos, last_pos)
    
    def overlap_in_first_pos(self, in_first_pos) :
        iso = self.input_selection_obj(in_first_pos)
        return iso.overlap_pos(in_first_pos)
        
    def overlap_in_first_last_pos(self, in_first_pos, in_last_pos) :
        iso = self.input_selection_obj(in_first_pos, in_last_pos)
        return iso.overlap_interval((in_first_pos, in_last_pos))
        
    def within_in_first_last_pos(self, in_first_pos, in_last_pos) :
        iso = self.input_selection_obj(in_first_pos, in_last_pos)
        return iso.within_in_first_last_pos(in_first_pos, in_last_pos)
        
    
    
    
    
    
    
    
    