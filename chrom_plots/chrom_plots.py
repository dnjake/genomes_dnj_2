# -*- coding: utf-8 -*-

from bokeh.plotting import figure
from bokeh.models import Range1d
from bokeh.layouts import column
from ..interval_series_plots import interval_gene_finder as igf

from . import chrom_event_data as ced


class gene_plot_cls(object) :
    plot_height = 100
    gene_y_bottom = 15
    gene_y_top = 65
    plot_y_top = 80
    selection_color = 'blue'
    def __init__(self, chrom) :
        self.chrom = chrom
        
    def do_plot(self, first, last, plot_width=1700, selection_range=None) :        
        x_axis_plot_width = plot_width
        x_pos_range = first, last
        genes_obj = igf.gene_finder_cls(self.chrom, x_pos_range, x_axis_plot_width,
                                                 self.gene_y_bottom, self.gene_y_top)
        genes_left = genes_obj.genes_left
        genes_top = genes_obj.genes_top
        genes_right = genes_obj.genes_right
        genes_bottom = genes_obj.gene_y_bottom
        gene_labels = genes_obj.gene_labels
        plot = figure(plot_width=x_axis_plot_width, plot_height=self.plot_height, 
                      tools=[], toolbar_location=None, title='genes')
        plot.yaxis.visible = None
        plot.yaxis.major_tick_line_color = None
        plot.yaxis.minor_tick_line_color = None
        plot.ygrid.grid_line_color = None
        plot.yaxis.major_label_text_font_size = '0pt'
        plot.toolbar.active_drag = None
        plot.y_range = Range1d(0, self.plot_y_top)  
        plot.x_range = Range1d(*x_pos_range)
        if genes_left is not None:
            plot.quad(left=genes_left, right=genes_right, top=genes_top,
                      bottom=genes_bottom, fill_color='white', line_color='grey')
        if selection_range is not None :
            selection_left, selection_right = selection_range
            left = [selection_left, selection_left]
            right = [selection_right, selection_right]
            top = [self.gene_y_bottom, self.plot_height]
            bottom = [0, self.gene_y_top]
            plot.quad(left=left, right=right, top=top, bottom=bottom, fill_color=self.selection_color,
                      line_color='white', line_width=0)
        if gene_labels is not None :
            plot.add_layout(gene_labels)
        return plot



class chrom_plot_cls(object) :
    def __init__(self, chrom) :
        self.chrom = chrom
        self.sro = ced.data_rdr_cls(chrom)
        self.gene_plot_obj = None
        
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
        p = figure(plot_width=800, plot_height=200, toolbar_location=None, title=field_name)
        p.y_range = Range1d(0, 1.1*float(y_max))
        p.x_range = Range1d(x_pos_first, x_pos_last)
        p.line(x, y)
        return p

    def do_column_plots(self, field_names, first=None, last=None, plot_width=800) :
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
            p = figure(plot_width=plot_width, plot_height=200, tools=[], toolbar_location=None, title=field_name)
            p.toolbar.active_drag = None
            p.y_range = Range1d(0, 1.1*float(y_max))
            p.x_range = Range1d(x_pos_first, x_pos_last)
            p.line(x, y)
            plots.append(p)
        p = plots[-1]
        p.xaxis.axis_label = 'Chromosome Position'
        return plots
            

    def do_standard_field_plots(self, first=None, last=None, plot_width=900) :
        field_names = ['samples_in_series', 'series_count', 'sample_weighted_length',
                       'sample_weighted_snps', 'snps_in_series', 'mean_series_length',
                       'mean_series_snps']
        return column(self.do_column_plots(field_names, first, last, plot_width))

    def do_chrom_plots(self, plot_width=900) :
        field_names = ['series_count', 'mean_series_length', 'sample_weighted_length', 'sample_weighted_snps']
        return column(self.do_column_plots(field_names, plot_width=plot_width))
    
    def do_simple_with_genes(self, first, last, plot_width=1700) :
        field_names = [ 'samples_in_series', 'series_count', 'mean_series_length']
        field_plots = self.do_column_plots(field_names, first, last, plot_width)
        if self.gene_plot_obj is None :
            self.gene_plot_obj = gene_plot_cls(self.chrom)
        gene_plot = self.gene_plot_obj.do_plot(first, last, plot_width)
        plots = [gene_plot]
        plots.extend(field_plots)
        return column(plots)
        
        

