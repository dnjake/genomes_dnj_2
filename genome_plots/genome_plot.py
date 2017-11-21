# -*- coding: utf-8 -*-

from . import genome_plot_data as gpd
from bokeh.plotting import figure
from bokeh.models import Range1d

        
def do_plot(field_name) :
    rdr = gpd.data_rdr_cls(field_name)
    d = rdr.genome_fraction_data()
    x = d['x']
    y = d['y']
    x_pos_first = x[0]
    x_pos_last = x[-1]
    y_max = y.max()
    p = figure(plot_width=800, plot_height=200, toolbar_location=None, title=field_name)
    p.y_range = Range1d(0, 1.1*float(y_max))
    p.x_range = Range1d(x_pos_first, x_pos_last)
    p.xaxis.axis_label = 'Fraction of Genome'
    p.line(x, y)
    return p
    
def do_log_plot(field_name) :
    rdr = gpd.data_rdr_cls(field_name)
    d = rdr.genome_fraction_data()
    x = d['x']
    y = d['y']
    m = y > 0
    x = x[m]
    y = y[m]
    #x_pos_first = x[0]
    x_pos_last = x[-1]
    #y_max = y.max()
    p = figure(plot_width=800, plot_height=200, toolbar_location=None, y_axis_type='log', title=field_name)
    #p.y_range = Range1d(0, 1.1*float(y_max))
    p.x_range = Range1d(0, x_pos_last)
    p.xaxis.axis_label = 'Fraction of Genome'
    p.line(x, y)
    return p
