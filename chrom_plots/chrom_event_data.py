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


stats_file_name = 'chrom_event_data.h5'
stats_table_name_end = '_event_data'
stats_file_path = os.path.join(mod_dir, stats_file_name)

def stats_table_name(chrom) :
    return 'chrom_' + str(chrom) + stats_table_name_end

class data_rdr_cls(object) :
    #x_axis_graph_count = 300000
    def __init__(self, chrom, x_samples=5000) :
        h5 = tb.open_file(stats_file_path, 'r')
        stats_table = getattr(h5.root, stats_table_name(chrom))
        self.data = stats_table[:]
        h5.close()
        self.chrom = chrom
        self.x_samples = x_samples
        self.event_last_pos = self.data['event_last_pos']
        self.chrom_first_pos = self.data['event_first_pos'][0]
        self.chrom_last_pos = self.data['event_last_pos'][-1]

    def x_interval(self, first_pos=None, last_pos=None) :
        if first_pos is None :
            first = self.chrom_first_pos
        else :
            first = first_pos
        if last_pos is None :
            last = self.chrom_last_pos
        else :
            last = last_pos
        first_f = float(first)            
        last_f = float(last)            
        interval = (last_f - first_f)/float(self.x_samples)
        sample_pos = np.arange(first_f, last_f, interval, dtype='f8')
        sample_pos += 0.5
        sample_pos = sample_pos.astype('u4')
        sample_pos[-1] = last
        sample_indexes = self.event_last_pos.searchsorted(sample_pos)
        return sample_indexes, sample_pos

    def data_dict_for_field_names(self, field_names, first_pos=None, last_pos=None) :
        sample_indexes, x_values = self.x_interval(first_pos, last_pos)
        out_dict = {'x': x_values}
        for name in field_names :
            out_dict[name] = self.data[name][sample_indexes]
        return out_dict


class chrom_event_data_writer_cls(object) :
    def __init__(self) :
        h5 = tb.open_file(fls.stats_file_path, 'r')
        stats_table = getattr(h5.root, fls.stats_table_name)
        self.event_data = stats_table[:]
        h5.close()

    def event_data_for_chrom(self, chrom) :
        bound = chrom + 1
        inds = self.event_data['chrom'].searchsorted([chrom, bound])
        event_data = self.event_data[inds[0]:inds[1]]
        return event_data

    def write_chrom_tables(self) :
        h5 = tb.open_file(stats_file_path, 'w', filters=h5_filters())
        for chrom in range(1, 23) :
            print 'processing', chrom
            chrom_data = self.event_data_for_chrom(chrom)
            stats_table = h5.create_table('/', stats_table_name(chrom), description=chrom_data.dtype)
            stats_table.append(chrom_data)
            stats_table.close()
        h5.close()
