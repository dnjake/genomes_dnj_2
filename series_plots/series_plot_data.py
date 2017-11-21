# -*- coding: utf-8 -*-

import tables as tb
import numpy as np
from ..stats_by_series import genome_by_series_stats as bss
from ..html_display import array_table as html



class series_plot_data_cls(object) :
    pops = ('no_afr', 'low_afr', 'some_afr', 'high_afr', 'all_afr')
    pop_types = [('no_afr', 'u4'), ('low_afr', 'u4'), ('some_afr', 'u4'), ('high_afr', 'u4'), ('all_afr', 'u4')]
    work_dtype = np.dtype(pop_types)
    x_axis_graph_count = 1000
    def __init__(self) :
        h5 = tb.open_file(bss.stats_file_path, 'r')
        stats_table = getattr(h5.root, bss.stats_table_name)
        self.data = stats_table[:]
        h5.close()
        top = float(self.data.size)
        interval = top/float(self.x_axis_graph_count)
        sample_indexes = np.arange(0.0, top, interval, dtype='f4')
        sample_indexes += 0.5
        self.sample_indexes = sample_indexes.astype('u4')
        self.sample_indexes[-1] = self.data.size - 1
        
    def find_pop_data(self, pt, pop_data) :
        m = self.data['pop_type'] == pt
        m = m.astype('u1')
        np.cumsum(m, dtype='u4', out=pop_data)

    def data_by_field_name(self, field_name) :
        self.data.sort(order=field_name)
        pop_type_series_counts = np.zeros(self.data.size, self.work_dtype)
        for pop in self.pops :            
            pt = bss.pop_types[pop]
            pop_data = pop_type_series_counts[pop]
            self.find_pop_data(pt, pop_data)
        data_type = []
        data_type.extend(self.pop_types)
        field_dtype = self.data[field_name].dtype
        field_type = ('field_value', field_dtype)
        data_type.extend([('all', 'u4'), field_type])
        data_dtype = np.dtype(data_type)
        self.field_data = np.zeros(self.sample_indexes.size, data_dtype)
        for pop in self.pops :
            self.field_data[pop] = pop_type_series_counts[pop][self.sample_indexes]
        self.field_data['field_value'] = self.data[field_name][self.sample_indexes]
        self.field_data['all'] = self.sample_indexes + 1
        return self.field_data
        










    '''
    this case is simpler because all I want is the offset indexes for the
    array.  I just need to use arange and the interval to find the index
    values into the array.  I need to turn those index values into fractions
    of the series for y coordinates and the series length for x coordinate
    
    This approach samples regularly by series number which may be the best
    The alternative is to sample by intervals of value
    
    What I want is a plot that shows the fraction of series less than a
    sample value.  I want to divide that fraction by the sub fraction
    of those series that belongs to different population types.  I want
    a graph that incrmentally adds the sub fraction for a pop type
    It would be nice if the segment of the graph for each pop type
    could be given a different background color.
    
    One way of doing the plot is drawing a patch for each pop type
    
    Another way is drawing a quad per pop type for each sample value
    probably could be a bar
    
    probably the best way is to generate the right rgb values for the
    plot
    
    brewer.py in the gallery looks like it is going what I want with pathces
    '''

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
