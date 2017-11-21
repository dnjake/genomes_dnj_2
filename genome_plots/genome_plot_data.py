# -*- coding: utf-8 -*-

import numpy as np
import tables as tb
import os

from ..stats_by_pos import genome_by_series_first_last_stats as fls

s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

def h5_filters() :
    return tb.Filters(complib='zlib', complevel=5)


stats_file_name = 'genome_sorted_stats_data.h5'
stats_table_name_end = '_sorted_stats'
stats_file_path = os.path.join(mod_dir, stats_file_name)

def stats_table_name(stat_name) :
    return stat_name + stats_table_name_end
    
class data_rdr_cls(object) :
    x_axis_graph_count = 1000
    def __init__(self, field_name) :
        h5 = tb.open_file(stats_file_path, 'r')
        data_table = getattr(h5.root, stats_table_name(field_name))
        self.data = data_table[:]
        h5.close()

    def genome_plot_data(self):
        data = self.data
        top_data = data[-1]
        top_offset = top_data['offset'] + top_data['length']
        start = data['offset'][0]
        interval = float(top_offset)/float(self.x_axis_graph_count)
        top_offset += (interval/2)
        sample_offsets = np.arange(start, top_offset, interval, dtype='f4')
        sample_offsets = sample_offsets.astype('u4')
        sample_indexes = data['offset'].searchsorted(sample_offsets)
        sample_indexes[-1] = data.size - 1
        graph_data_type = np.dtype([('x', 'u4'), ('y', data['value'].dtype)])
        graph_data = np.zeros(sample_indexes.size, graph_data_type)
        graph_data['x'] = data['offset'][sample_indexes]
        graph_data['y'] = data['value'][sample_indexes]
        return graph_data

    def genome_fraction_data(self) :
        field_data = self.genome_plot_data()
        x = field_data['x']
        x_max = x.max()
        y = field_data['y']
        data_dtype = np.dtype([('x', 'f4'), ('y', y.dtype)])
        out_data = np.zeros(x.size, data_dtype)
        out_data['x'] = x.astype('f4')/float(x_max)
        out_data['y'] = y
        return out_data


class genome_field_data_writer(object) :
    centrosome_subtract_count = 23
    stats_fields = ('mean_series_length', 'samples_in_series', 'series_count', 'sample_weighted_length',
                    'snps_in_series', 'mean_series_snps', 'sample_weighted_snps', 'no_afr_series_count',
                    'no_afr_sample_weighted_length', 'low_afr_series_count', 'low_afr_sample_weighted_length',
                    'some_afr_series_count', 'some_afr_sample_weighted_length', 'high_afr_series_count',
                    'high_afr_sample_weighted_length', 'all_afr_series_count', 'all_afr_sample_weighted_length')
    def __init__(self) :
        h5 = tb.open_file(fls.stats_file_path, 'r')
        stats_table = getattr(h5.root, fls.stats_table_name)
        self.event_data = stats_table[:]
        h5.close()
        self.subtract_centrosomes()
        
    def subtract_centrosomes(self) :
        d = self.event_data
        lengths = d['event_last_pos'] - d['event_first_pos']
        sa = lengths.argsort()
        m = np.ones(d.size,'?')
        m[sa[-self.centrosome_subtract_count:]] = False
        self.event_data = self.event_data[m]
        self.event_data_lengths = self.event_data['event_last_pos'] - self.event_data['event_first_pos']
        
        
    def generate_offsets(self, data) :
        np.cumsum(data['length'], dtype='u4', out=data['offset'])
        '''
        offset = 0 
        for item in data :
            item['offset'] = offset
            offset += item['length']
        '''

    def length_value_dtype(self, field_data) :
        return np.dtype([('offset', 'u4'), ('length', 'u4'), ('value', field_data.dtype)])
        
    def generate_value_data(self, field_name) :
        field_data = self.event_data[field_name]
        value_data = np.zeros(field_data.size, self.length_value_dtype(field_data))
        value_data['value'] = field_data
        value_data['length'] = self.event_data['event_last_pos'] - self.event_data['event_first_pos']
        return value_data
    
    def genome_data_for_field(self, field_name) :
        value_data = self.generate_value_data(field_name)
        value_data.sort(order='value')
        self.generate_offsets(value_data)
        return value_data

    def write_field_tables(self) :
        h5 = tb.open_file(stats_file_path, 'w', filters=h5_filters())
        for field_name in self.stats_fields :
            print 'processing', field_name
            field_data = self.genome_data_for_field(field_name)
            stats_table = h5.create_table('/', stats_table_name(field_name), description=field_data.dtype)
            stats_table.append(field_data)
            stats_table.close()
        h5.close()

