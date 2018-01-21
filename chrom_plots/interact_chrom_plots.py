# -*- coding: utf-8 -*-

from bokeh.layouts import layout, column
from bokeh.models import ColumnDataSource, Range1d, HoverTool
import genomes_dnj_2.chrom_plots.chrom_event_data as ced

from bokeh.plotting import figure


class chrom_interact_cls(object) :
    hover = HoverTool(tooltips = [
                            ('index', '$index'),
                            ('x', '$x'),
                            ('pos', '@x')
                            ])
                            
    tools_x = 'pan, xzoom_in, xzoom_out'
    
    tools_y = 'pan, yzoom_in, yzoom_out'

    data_items = ['samples_in_series', 'series_count', 'mean_series_length'] 
    x_samples = 300000
    def __init__(self, chrom, plot_width=1600) :
        self.chrom = chrom
        self.plot_width = plot_width
        self.sro = ced.data_rdr_cls(self.chrom, self.x_samples)
        
    def do_plots(self, first=None, last=None) :
        data_items = self.data_items
        source_dict = self.sro.data_dict_for_field_names(data_items)
        source = ColumnDataSource(source_dict)        
        plots = []
        for ind in range(len(data_items)) :
            if ind == 0 :
                tools = self.tools_x
            else :
                tools = self.tools_y
            f = data_items[ind]
            p = figure(plot_width=self.plot_width, plot_height=250, title=f, tools=tools)
            p.line('x', f, source=source)
            plots.append(p)
        if first is None :
            x_min = 0
        else :
            x_min = first
        if last is None :
            x_max = source.data['x'][-1]
        else :
            x_max = last
        x_range = Range1d(x_min, x_max)
        for p in plots[:] :
            p.x_range = x_range
        p0 = plots[0]
        p0.add_tools(self.hover)
        p0.toolbar.active_inspect = None
        return column(plots)