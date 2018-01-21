# -*- coding: utf-8 -*-

import numpy as np
from bokeh.plotting import figure
from bokeh.models import Range1d

class snp_match_count_plot_cls(object) :
    line_data_dtype = np.dtype([('x_snp_count', 'f4'), ('y_chrom_count', 'f4')])
    plot_height = 400
    def __init__(self, count_data, plot_width=900, toolbar_location=None) :
        self.line_coords_from_unique_counts(count_data)
        self.plot_width = plot_width
        self.toolbar_location = toolbar_location

    def line_coords_from_unique_counts(self, count_data) :
        line_coords = [(-0.5,0)]
        x_last = -0.5
        for val, count in count_data :
            val -= 0.5
            if val == x_last :
                line_coords.append((x_last, count))
                x_last += 1
                line_coords.append((x_last, count))
            else :
                line_coords.append((x_last, 0))
                x_last = val
                line_coords.append((x_last, 0))
                line_coords.append((x_last, count))
                x_last += 1
                line_coords.append((x_last, count))
        line_coords.append((x_last, 0))
        #print line_coords
        self.x_y_data = np.array(line_coords, self.line_data_dtype)
                
    def do_plot(self) :
        self.plot_figure = figure(plot_width=self.plot_width, plot_height=self.plot_height,
                                  toolbar_location=self.toolbar_location)
        x_vals = self.x_y_data['x_snp_count']
        y_vals = self.x_y_data['y_chrom_count']
        x_max = x_vals.max() + 1
        y_max = 1.1*y_vals.max()
        self.plot_figure.y_range = Range1d(0, y_max)
        self.plot_figure.x_range = Range1d(-0.5, x_max)
        self.plot_figure.line(x=x_vals, y=y_vals)
        return self.plot_figure
        
        
        
        
