# -*- coding: utf-8 -*-

import numpy as np
import tables as tb
from . import genome_plot_data as gpd
from ..html_display import array_table as html


'''
want to get means, medians, and max for each of the key statistics
mean is just the value that divides the genome in half
also want a bin table for each statistic that shows the fraction
of the genome with a larger or equal value to the bin

want to move subtracting the centrosomes and probably opening
the table to an helper module.  May also want to put functions
for converting values into that module

probably also could vectorize the code for generating plot values
should be able to use np.cumsum
'''

names_ordered_key_stats = ( 'samples_in_series', 'series_count', 'sample_weighted_length', 'sample_weighted_snps',
                            'snps_in_series',   'mean_series_length', 'mean_series_snps')

samples_in_series_bins = np.array([1000, 2000, 3000, 4000, 4500], dtype='u2')
series_count_bins = np.array([5, 10, 20, 40, 80],dtype='u2')
sample_weighted_length_bins = np.array([5e7, 1e8, 2e8, 4e8, 8e8, 1.6e9, 3.2e9, 6.4e9], dtype='u4')
sample_weighted_snps_bins = np.array([5e3, 1e4, 2e4, 4e4, 8e4, 1.6e5, 3.2e5], dtype='u4')
snps_in_series_bins = np.array([50, 100, 200, 400, 800, 1600], dtype='u2')
mean_series_length_bins = np.array([10000, 20000, 40000, 80000, 160000, 320000, 640000], dtype='u4')
mean_series_snps_bins = np.array([5, 10, 15, 20, 25, 30], dtype='u2')

field_bins = { 'samples_in_series': samples_in_series_bins,
               'series_count': series_count_bins,
               'sample_weighted_length': sample_weighted_length_bins,
               'sample_weighted_snps': sample_weighted_snps_bins,
               'snps_in_series': snps_in_series_bins,
               'mean_series_length' : mean_series_length_bins,
               'mean_series_snps' : mean_series_snps_bins
            }        
class field_event_table_stats_cls(object) :
    def __init__(self, field_name) :
        self.field_name = field_name
        h5 = tb.open_file(gpd.stats_file_path, 'r')
        data_table = getattr(h5.root, gpd.stats_table_name(field_name))
        self.data = data_table[:]
        h5.close()
        self.offset_data = self.data['offset']
        self.length_data = self.data['length']
        self.value_data = self.data['value']
        
    def field_mean(self) :
        values = self.value_data.astype('f4')
        lengths = self.length_data.astype('f4')
        weighted = values*lengths
        mean = weighted.sum()/lengths.sum()
        return mean
        
    def field_median(self) :
        offset_max = self.offset_data[-1] + self.length_data[-1]
        offset_median = offset_max/2
        ind_median = self.offset_data.searchsorted(offset_median)
        offset_high = self.offset_data[ind_median] 
        value_high = self.value_data[ind_median]
        if offset_high == offset_median :
            return value_high 
        else :
            ind_low = ind_median - 1
            offset_low = self.offset_data[ind_median]
            value_low = self.value_data[ind_low]
            delta = (value_high-value_low)*((offset_median-offset_low)/(offset_high-offset_low))
            return value_low + delta        
         
    def field_fraction_bins(self) :
        bins = field_bins[self.field_name]
        inds = self.value_data.searchsorted(bins)
        offset_values = self.offset_data[inds]
        offset_max = float(self.offset_data[-1])
        offset_values = offset_values.astype('f4')
        offset_values = offset_max - offset_values
        offset_values = offset_values/offset_max
        bin_data_dtype = np.dtype([('data_value', bins.dtype), ('genome_fract', 'f4')])
        bin_data = np.zeros(offset_values.size, bin_data_dtype)
        bin_data['data_value'] = bins
        bin_data['genome_fract'] = offset_values
        return bin_data
            
        
class genome_stats_html_cls(object) :
        
    def html_from_field_data(self, field_name, field_mean, field_median, field_data) :
        fmt = html.big_fmt
        if field_name == 'sample_weighted_length' :
            fmt = html.exp_fmt
        h = ['<tr><th>mean</th><th>median</th>']
        r = ['<tr>']
        r.append(html.num_tag + fmt.format(int(field_mean)) + '</td>')
        r.append(html.num_tag + fmt.format(int(field_median)) + '</td>')
        for hv, dv in field_data :
            h.append('<th>')
            h.append(fmt.format(hv))
            h.append('</th>')
            r.append(html.num_tag)
            r.append(html.float_fmt.format(dv))
            r.append('</td>')
        h.append('</tr>')
        r.append('</tr>')
        html_table = ['<h4>' + field_name + '</h4>']
        html_table.append('<table>')
        html_table.append(''.join(h))
        html_table.append(''.join(r))
        html_table.append('</table>')
        table_html = '\n'.join(html_table)
        return table_html       
        
    def html_for_field(self, field_name) :
        fso = field_event_table_stats_cls(field_name)
        mean = fso.field_mean()
        median = fso.field_median()
        data = fso.field_fraction_bins()
        return self.html_from_field_data(field_name, mean, median, data)
        
        
    def genome_stats_html(self) :
        stats_html = []
        for field_name in names_ordered_key_stats :
            stats_html.append(self.html_for_field(field_name))
        return '\n'.join(stats_html)
        
        
        
        
        
        
        
         
         
         
         