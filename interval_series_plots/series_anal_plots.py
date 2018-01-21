# -*- coding: utf-8 -*-

#import numpy as np
from bokeh.layouts import column
from bokeh.models import Div
from . import plot_style as style
from ..chrom_plots import chrom_plots as cp
from ..interval import series_select as ssel
from . import series_sample_plot_layout as spl
from . import series_sample_data_table as sdt
import genomes_dnj_2.interval.series_input as sin

def add_to_layout(layout_items, series_plot_and_data_table) :
    plot, table = series_plot_and_data_table
    layout_items.append([plot])
    layout_items.append([Div(text=table)])

class series_anal_plt_cls(object):
    series_stats_names = [ 'samples_in_series', 'series_count', 'mean_series_length']
    plot_style = style
    def __init__(self, series_anal_obj) :
        self.chrom = series_anal_obj.chrom        
        self.series_anal_obj = series_anal_obj
        self.anal_pos_range = self.series_anal_obj.anal_first_pos, self.series_anal_obj.anal_last_pos
        self.anal_selection_range = self.series_anal_obj.selection_range
        self.selection_obj = ssel.select_series_cls(self.series_anal_obj)
        first_pos, last_pos = self.anal_selection_range
        self.interval_id = str(first_pos) + '_' + str(last_pos)
        
    @classmethod
    def all_read_series_from_chrom_data(cls, chrom, first_pos, last_pos) :
        in_rdr = sin.series_data_rdr_cls(chrom)
        in_data_obj = in_rdr.input_selection_obj(first_pos, last_pos)
        min_first_pos = in_data_obj.series_first_pos.min()
        max_last_pos = in_data_obj.series_last_pos.max()
        series_anal_obj = in_data_obj.overlap_interval((min_first_pos, max_last_pos))
        return cls(series_anal_obj)
    
    @classmethod 
    def selected_series_from_chrom_data(cls, chrom, first_pos, last_pos) :
        in_rdr = sin.series_data_rdr_cls(chrom)
        in_data_obj = in_rdr.input_selection_obj(first_pos, last_pos)
        series_anal_obj = in_data_obj.overlap_interval((first_pos, last_pos))
        return cls(series_anal_obj)
    
    @classmethod
    def selected_series_from_input_data(cls, in_data_obj, first_pos, last_pos) :
        series_anal_obj = in_data_obj.overlap_interval((first_pos, last_pos))
        return cls(series_anal_obj)

    def do_series_plot(self, po) :
        po.plot_data = self.series_anal_obj.plot_data_from_series_sample_match_data(po.series_sample_match_data)
        has_matches = True
        if po.plot_data['match_allele_count'] is None :
            has_matches = False
        height_from_matches = True
        if po.plot_series_height_from_series_allele_count :
            height_from_matches = False
        plt_obj = spl.series_plot_cls(self.chrom, self.anal_pos_range, po.plot_data,
                                      x_axis_screen_width=self.plot_style.plot_width,
                                      title=po.plot_title, has_matches=has_matches,
                                      height_from_matches=height_from_matches)
        plt_obj.do_plot()
        po.series_plot = plt_obj.plot

    def do_series_data_table(self, po) :
        sdt_obj = sdt.series_data_table_cls(po)
        po.series_data_table_html = sdt_obj.series_data_html()

    def do_common(self, po) :
        self.do_series_plot(po)
        if po.do_gene_plot :
            self.do_interval_gene_plot(po)
        if po.do_stats_plots :
            self.do_interval_stats_plots(po)
        if po.do_series_data_table :
            self.do_series_data_table(po)
        po.assemble_layout()

    def get_interval_gene_plot(self, interval_obj, selection_range=None) :
        gene_plot_obj = cp.gene_plot_cls(interval_obj.chrom)
        first, last = interval_obj.anal_pos_range
        plot_width = self.plot_style.plot_width
        gene_plot = gene_plot_obj.do_plot(first, last, plot_width, selection_range=selection_range)
        return gene_plot    
        
    def get_interval_series_stats_plots(self, interval_obj, stats_names) :
        chrom_plot_obj = cp.chrom_plot_cls(interval_obj.chrom)        
        first, last = interval_obj.anal_pos_range
        stats_plots = chrom_plot_obj.do_column_plots(stats_names, first, last, self.plot_style.plot_width)
        return stats_plots
        
        
    def get_interval_gene_and_series_stats_plots(self, interval_obj) :
        plots = [self.get_interval_gene_plot(interval_obj)]
        plots.extend(self.get_interval_series_stats_plots(interval_obj))
        return plots

    def do_interval_gene_plot(self, po) :
        if (po.do_selection_range is None) or (not po.do_selection_range) :
            po.gene_plot = self.get_interval_gene_plot(self)
        else :
            po.gene_plot = self.get_interval_gene_plot(self, selection_range=self.anal_selection_range)
    
    def do_interval_stats_plots(self, po) :
        po.stats_plots = self.get_interval_series_stats_plots(self, po.stats_plot_names)
        
    def do_all_series_plot(self, po) :
        po.plot_title = 'all_series_' + self.interval_id
        po.series_sample_match_data = self.series_anal_obj.out_series_sample_data()
        self.do_common(po)
        return po

    
    def do_most_common_series_plot(self, po) :
        po.series_sample_match_data = self.selection_obj.select_by_most_common()
        po.plot_title = 'most_common_series_' + self.interval_id
        self.do_common(po)
        return po
    
    def do_basis_series_plot(self, po) :
        po.select_series_data_index = None
        po.select_allele_mask = None
        if po.select_series_index_and_id is not None :
            po.select_series_data_index, series_id = po.select_series_index_and_id
            po.plot_title = 'basis_series_' + series_id
        elif po. select_allele_mask_and_id is not None  :
            po.select_allele_mask, allele_mask_id = po.select_allele_mask_and_id
            po.plot_title = 'basis_series_' + allele_mask_id
        else :
            po.plot_title = 'basis_series_' + self.interval_id
        self.selection_obj.select_by_basis_from_plt_obj(po)
        self.do_common(po)
        return po

    def do_longest_series_plot(self, po) :
        po.series_sample_match_data = self.selection_obj.select_by_length()
        po.plot_title = 'longest_series_' + self.interval_id
        self.do_common(po)
        return po

    def add_min_match_to_plot_title(self, po) :
        min_match_str = self.plot_style.float_fmt.format(po.min_match)
        po.plot_title = po.plot_title + '_min_match_' + min_match_str
        
    def do_hierarchy_series_plot(self, po) :
        series_di, series_id = po.hierarchy_series_index_and_id
        po.hierarchy_data_index = series_di
        self.series_anal_obj.hierarchy_data_from_plt_obj(po)
        po.plot_title = 'hierarchy_' + series_id + '_data_index_' + str(series_di)
        self.add_min_match_to_plot_title(po)
        self.do_common(po)
        return po
        
    def do_yes_no_id(self, po) :
        if po.yes_series_data_indexes_and_ids is None :
            po.yes_series_data_indexes = None
            po.yes_series_ids = None
            po.yes_no_id = None            
        else :
            po.yes_series_data_indexes, po.yes_series_ids = zip(*po.yes_series_data_indexes_and_ids)
            sids = '_'.join(po.yes_series_ids)
            po.yes_no_id = 'yes_' + sids
        if po.no_series_data_indexes_and_ids is None :
            po.no_series_data_indexes = None
            po.no_series_ids = None
        else :
            po.no_series_data_indexes, po.no_series_ids = zip(*po.no_series_data_indexes_and_ids)            
            nids = '_'.join(po.no_series_ids)
            if po.yes_no_id is None :
                po.yes_no_id = 'no_' + nids
            else :
                po.yes_no_id = po.yes_no_id + '_no_' + nids

    def do_yes_no_series_plot(self, po) :
        self.do_yes_no_id(po)
        po.plot_title = po.yes_no_id
        self.add_min_match_to_plot_title(po)
        self.series_anal_obj.sub_super_yes_no_from_plt_obj(po)
        self.do_common(po)
        return po

    def do_allele_mask_yes_no_series_plot(self, po) :
        self.do_yes_no_id(po)
        po.yes_allele_mask, po.allele_mask_id = po.yes_allele_mask_and_id
        if po.yes_no_id is None :
            po.plot_title = po.allele_mask_id
        else :
            po.plot_title = po.allele_mask_id + '_' + po.yes_no_id
        self.add_min_match_to_plot_title(po)
        self.series_anal_obj.sub_super_allele_mask_yes_no_from_plt_obj(po)
        self.do_common(po)
        return po
    
    











        