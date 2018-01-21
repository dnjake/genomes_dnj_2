# -*- coding: utf-8 -*-
from bokeh.models import Range1d
from bokeh.plotting import figure

from . import series_row_layout as srl
from . import series_color_finder as scf
from . import series_snp_finder as ssnpf
from . import series_height_finder as shf
from . import series_y_values_finder as yvf
from . import series_id_labels as sil

class series_plot_cls(object) :
    
    def __init__(self, chrom, x_pos_range, series_plot_data, x_axis_screen_width=900, title=None,
                 has_matches=True, height_from_matches=False, color_from_matches=False) :
        self.chrom = chrom
        self.x_pos_range = x_pos_range        
        self.series_plot_data = series_plot_data
        self.x_axis_screen_width = x_axis_screen_width
        self.title = title
        self.has_matches = has_matches
        self.height_from_matches = height_from_matches
        self.color_from_matches = color_from_matches
        self.data_index = series_plot_data['data_index']
        self.first_pos = series_plot_data['first_pos']
        self.last_pos = series_plot_data['last_pos']
        self.snp_count = series_plot_data['snp_count']
        self.allele_count = series_plot_data['allele_count']
        self.allele_mask = series_plot_data['allele_mask']
        self.series_data = series_plot_data['series_data']
        if self.has_matches :
            self.match_allele_count = series_plot_data['match_allele_count']
            self.match_allele_mask = series_plot_data['match_allele_mask']
        self.height_allele_count = self.allele_count
        self.order_allele_count = self.allele_count
        if self.height_from_matches :
            self.height_allele_count = self.match_allele_count
        
    def do_row_layout(self) :
        self.row_layout_obj = srl.row_layout_cls(self.first_pos, self.last_pos, self.height_allele_count, self.order_allele_count,
                                                 self.x_pos_range, self.x_axis_screen_width, self.has_matches)
        self.layout_rows = self.row_layout_obj.layout_rows
        self.row_height_allele_count = self.layout_rows['height_allele_count']
        self.series_row_layout = self.row_layout_obj.series_row_layout
        self.series_row_index = self.series_row_layout['row_index']
        self.series_space_start = self.series_row_layout['pos_start']
        self.series_space_end = self.series_row_layout['pos_end']

    def do_series_colors(self) :
        allele_mask = self.allele_mask
        if self.color_from_matches :
            allele_mask = self.match_allele_mask
        self.series_colors_obj = scf.color_finder_cls(allele_mask)
        self.series_color = self.series_colors_obj.series_color

    def do_series_snps(self) :
        self.series_snp_finder_obj = ssnpf.snp_finder_cls(self.chrom, self.series_data)
        self.series_snps_pos = self.series_snp_finder_obj.series_snps_pos
        
    def do_series_heights(self) :
        self.height_finder_obj = shf.row_series_height_finder_cls(self.height_allele_count,
                                                                  self.row_height_allele_count)
        self.row_height = self.height_finder_obj.row_height
        self.series_height = self.height_finder_obj.series_height
        
    def do_series_y_values(self) :
        self.y_values_obj =  yvf.y_values_finder_cls(self.row_height,
                                        self.series_row_index, self.series_height)
        self.series_y_bottom = self.y_values_obj.series_y_values['y_bottom'].astype('i4')
        self.series_y_top = self.y_values_obj.series_y_values['y_top'].astype('i4')
        self.plot_top = int(self.y_values_obj.plot_top)

    def do_series_id_labels(self) :
        self.id_labels_obj = sil.id_labels_cls(self.snp_count, self.allele_count,
                                               self.series_space_start, self.series_y_top)
        self.series_id_labels = self.id_labels_obj.id_labels
        
         
    def do_snp_lines(self) :        
        self.snp_lines_obj = ssnpf.snp_line_plot_cls(self.series_snps_pos, self.series_y_bottom,
                                                     self.series_y_top)
        self.snp_x_vals = self.snp_lines_obj.snp_x_vals
        self.snp_y_vals = self.snp_lines_obj.snp_y_vals
        
    def do_series_plot(self) :
        self.plot = figure(plot_width=self.x_axis_screen_width, plot_height=self.plot_top,
                           tools=[], toolbar_location=None, title=self.title)
        self.plot.quad(left=self.first_pos, right=self.last_pos, top=self.series_y_top, 
                       bottom=self.series_y_bottom, fill_color=self.series_color, alpha=1.0, line_color='black' )
        self.plot.add_layout(self.series_id_labels)
        self.plot.multi_line(self.snp_x_vals, self.snp_y_vals, color='black')
        self.plot.yaxis.visible = None
        #self.plot.yaxis.axis_line_width = 0
        self.plot.yaxis.major_tick_line_color = None
        self.plot.yaxis.minor_tick_line_color = None
        self.plot.yaxis.major_label_text_font_size = '0pt'
        self.plot.toolbar.active_drag = None
        #self.plot.toolbar.active_drag = None
        self.plot.y_range = Range1d(0, self.plot_top)
        self.plot.x_range = Range1d(*self.x_pos_range)
                
    def do_series_match_labels(self) :
        self.match_labels_obj = sil.match_labels_cls(self.match_allele_count, self.series_space_end, self.series_y_top)
        self.series_match_labels = self.match_labels_obj.match_labels
        self.plot.add_layout(self.series_match_labels)


    def do_plot(self) :
        self.do_row_layout()
        self.do_series_colors()
        self.do_series_snps()
        self.do_series_heights()
        self.do_series_y_values()
        self.do_series_id_labels()
        self.do_snp_lines()
        self.do_series_plot()
        if self.has_matches :
            self.do_series_match_labels()
        
        
        
        
        
        
        
        
        
