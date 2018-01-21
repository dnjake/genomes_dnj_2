# -*- coding: utf-8 -*-
import tables as tb
import numpy as np
from ..gene_data import gene_data_rdr as gdr
from bokeh.layouts import layout, column
from bokeh.models import ColumnDataSource, Range1d, HoverTool, LabelSet
import genomes_dnj_2.chrom_plots.chrom_event_data as ced

from bokeh.plotting import figure


class chrom_genes_cls(object) :
    plot_height = 100
    gene_y_bottom = 15
    gene_y_top = 65
    gene_mid_y = (gene_y_bottom + gene_y_top)/2
    plot_y_top = 80
    min_data_space_per_gene_letter = 5000
    max_sym = 10
    max_bound = max_sym - 1
    sym_dtype = np.dtype([('mid_pos', 'u4'), ('sym', 'S10')])
    def __init__(self, chrom) :
        self.chrom = chrom
        gro = gdr.genes_rdr_cls(self.chrom)
        self.genes = gro.genes
        self.pos_left = self.genes['tx_start']
        self.pos_right = self.genes['tx_end']

    def find_gene_syms(self) :
        sym_x_coord = []
        last_end = 0
        for chrom, tx_start, tx_end, sym in self.genes :
            sym_len = len(sym)
            if sym_len > self.max_sym :
                sym = sym[:self.max_bound] + '-'
                sym_len = 10
            sym_mid = (tx_start + tx_end)/2
            half_sym_len = sym_len/2
            half_sym_space = self.min_data_space_per_gene_letter*half_sym_len
            if sym_mid - half_sym_space > last_end :
                sym_x_coord.append((sym_mid, sym))
                last_end = sym_mid + half_sym_space
        self.text_coord = np.array(sym_x_coord, dtype=self.sym_dtype)
        self.text_mid_pos = self.text_coord['mid_pos']
        self.text_sym = self.text_coord['sym']

    def do_plot(self, screen_width=1700) :
        self.find_gene_syms()
        x_axis_screen_width = screen_width
        genes_left = self.pos_left
        genes_top = self.gene_y_top
        genes_right = self.pos_right
        genes_bottom = self.gene_y_bottom
        tools_x = 'pan, xzoom_in, xzoom_out'
        plot = figure(plot_width=x_axis_screen_width, plot_height=self.plot_height, 
                      tools=tools_x, toolbar_location=None, title='genes')
        plot.yaxis.visible = None
        plot.yaxis.major_tick_line_color = None
        plot.yaxis.minor_tick_line_color = None
        plot.ygrid.grid_line_color = None
        plot.yaxis.major_label_text_font_size = '0pt'
        plot.y_range = Range1d(0, self.plot_y_top)  
        plot.quad(left=genes_left, right=genes_right, top=genes_top,
                  bottom=genes_bottom, fill_color='white', line_color='grey')
        mid_y = np.zeros(self.text_mid_pos.size,'i4')
        mid_y[:] = self.gene_mid_y
        self.data_source = ColumnDataSource({'x': self.text_mid_pos, 'y': mid_y, 'vals': self.text_sym})                
        self.gene_labels = LabelSet(x='x', y='y', text='vals', source=self.data_source, level='glyph',
                                    render_mode='canvas', text_baseline='middle', text_align='center',
                        text_font_size=('8pt'), text_font_style=('bold'), text_alpha=1.0)
        plot.add_layout(self.gene_labels)
        return plot


        
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
        self.genes = chrom_genes_cls(self.chrom)
        
    def do_plots(self) :
        data_items = self.data_items
        source_dict = self.sro.data_dict_for_field_names(data_items)
        source = ColumnDataSource(source_dict)        
        plots = []
        gene_plot = self.genes.do_plot(self.plot_width)
        plots.append(gene_plot)
        for ind in range(len(data_items)) :
            if ind == 0 :
                tools = self.tools_x
            else :
                tools = self.tools_y
            f = data_items[ind]
            p = figure(plot_width=self.plot_width, plot_height=250, title=f, tools=tools)
            p.line('x', f, source=source)
            plots.append(p)
        x_min = 0
        x_max = source.data['x'][-1]
        x_range = Range1d(x_min, x_max)
        for p in plots[:] :
            p.x_range = x_range
        p1 = plots[1]
        p1.add_tools(self.hover)
        p1.toolbar.active_inspect = None
        return column(plots)