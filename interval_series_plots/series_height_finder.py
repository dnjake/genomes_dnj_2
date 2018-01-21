# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np

class row_series_height_finder_cls(object) :
    allele_mult = 6.0                                
    def __init__(self, series_allele_count, row_allele_count) :
        self.series_allele_count = series_allele_count
        self.row_allele_count = row_allele_count
        self.series_height = np.zeros(self.series_allele_count.size, dtype='f4')
        self.row_height = np.zeros(self.row_allele_count.size)
        self.calc_heights()
        
    def calc_heights(self) :
        self.series_height = np.log2(self.series_allele_count) 
        mask = self.series_height >= 4.0
        self.series_height[mask] = self.series_height[mask] - 3.0
        mask = np.logical_not(mask)
        self.series_height[mask] = 1.0
        self.series_height = self.allele_mult*self.series_height
        self.row_height = np.log2(self.row_allele_count) 
        mask = self.row_height >= 4.0
        self.row_height[mask] = self.row_height[mask] - 3.0
        mask = np.logical_not(mask)
        self.row_height[mask] = 1.0
        self.row_height = self.allele_mult*self.row_height
        
        