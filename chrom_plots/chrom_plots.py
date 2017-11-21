# -*- coding: utf-8 -*-

from bokeh.plotting import figure
from bokeh.models import Range1d
from bokeh.layouts import column

from . import chrom_event_data as ced

class chrom_plot_cls(object) :
    def __init__(self, chrom) :
        self.sro = ced.data_rdr_cls(chrom)
        
    def do_plot(self, field_name, first=None, last=None) :
        plot_dict = self.sro.data_dict_for_field_names([field_name], first, last)
        x = plot_dict['x']
        y = plot_dict[field_name]
        if first is None :
            x_pos_first = 0
        else :
            x_pos_first = plot_dict['x'][0]
        x_pos_last = x[-1]
        y_max = y.max()
        p = figure(plot_width=900, plot_height=250, toolbar_location=None, title=field_name)
        p.y_range = Range1d(0, 1.1*float(y_max))
        p.x_range = Range1d(x_pos_first, x_pos_last)
        p.line(x, y)
        return p

    def do_column_plots(self, field_names, first=None, last=None, plot_width=900) :
        plot_dict = self.sro.data_dict_for_field_names(field_names, first, last)
        x = plot_dict['x']
        if first is None :
            x_pos_first = 0
        else :
            x_pos_first = plot_dict['x'][0]
        x_pos_last = x[-1]
        plots = []
        for field_name in field_names :
            y = plot_dict[field_name]
            y_max = y.max()
            p = figure(plot_width=plot_width, plot_height=250, toolbar_location=None, title=field_name)
            p.y_range = Range1d(0, 1.1*float(y_max))
            p.x_range = Range1d(x_pos_first, x_pos_last)
            p.line(x, y)
            plots.append(p)
        p = plots[-1]
        p.xaxis.axis_label = 'Chromosome Position'
        return column(plots)
            

    def do_standard_field_plots(self, first=None, last=None, plot_width=900) :
        field_names = ['samples_in_series', 'series_count', 'sample_weighted_length',
                       'sample_weighted_snps', 'snps_in_series', 'mean_series_length',
                       'mean_series_snps']
        return self.do_column_plots(field_names, first, last, plot_width)

    def do_chrom_plots(self, plot_width=900) :
        field_names = ['series_count', 'mean_series_length', 'sample_weighted_length', 'sample_weighted_snps']
        return self.do_column_plots(field_names, plot_width=plot_width)

