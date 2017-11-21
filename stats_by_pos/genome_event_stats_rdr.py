# -*- coding: utf-8 -*-

import numpy as np
import tables as tb
from . import genome_by_series_first_last_stats as fls
from . import by_pos_stats_helpers as psh

class stats_rdr_cls(object) :
    #pos_value_fields = [('event_first_pos', '<u4'), ('event_last_pos', '<u4')]
    #, ('value', 'u2')])
    #graph_data_fields = [('x', 'u4')]
    #, ('y', 'u2')])
    def __init__(self, chrom=None) :
        h5 = tb.open_file(fls.stats_file_path, 'r')
        stats_table = getattr(h5.root, fls.stats_table_name)
        self.event_data = stats_table[:]
        h5.close()
        self.chrom = chrom
        if self.chrom is not None :
            bound = chrom + 1
            inds = self.event_data['chrom'].searchsorted([self.chrom, bound])
            self.event_data = self.event_data[inds[0]:inds[1]]

    def plot_data_for_field(self, field_name) :
        return psh.graph_data_by_field(self.event_data, field_name)       

