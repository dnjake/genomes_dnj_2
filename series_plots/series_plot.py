# -*- coding: utf-8 -*-

import numpy as np
from . import series_plot_data as spd
from bokeh.plotting import figure
from bokeh.palettes import brewer
from bokeh.models import ColumnDataSource, Range1d

class data_plot_cls(object) :
    def __init__(self) :
        self.data_rdr = spd.series_plot_data_cls()
        
    def do_log_plot(self, field_name) :
        d = self.data_rdr.data_by_field_name(field_name)
        pops = self.data_rdr.pops
        pop_data = {}
        all_d = d['all']
        all_max = float(all_d[-1])
        for pop in pops :
            pd = d[pop].astype('f4')
            #pd_max = pd[-1]
            #pd = pd_max - pd
            pd = pd/all_max
            pop_data[pop] = pd
        areas = {}
        pd_last = np.zeros(all_d.size, dtype='f4')
        for pop in pops :
            pd = pop_data[pop]
            pd += pd_last
            areas[pop] = np.concatenate((pd_last[::-1], pd))
            pd_last = pd
        x = d['field_value']
        x_min = x[0]
        x_max = x[-1]
        x2 = np.concatenate((x[::-1], x))
    
        px = []
        py = []
        for pop in pops :
            px.append(x2)
            py.append(areas[pop])
        colors = brewer["Spectral"][len(pops)]
        source = ColumnDataSource(dict(
            xs=px,
            ys=py,
            color=colors,
            label=['no afr', 'low afr', 'some afr', 'high afr', 'all afr']
        ))
            
        x_axis_label = 'series ' + field_name
        p = figure(plot_width=900, plot_height=500, toolbar_location=None, title=field_name,
                   x_axis_type='log', x_range=Range1d(x_min, x_max))
                   
        #print colors
        p.patches( xs='xs', ys='ys', color='color', legend='label', source=source, line_color=None)
        p.legend.location = 'top_left'
        p.yaxis.axis_label = 'fraction of series'
        
        return p

