# -*- coding: utf-8 -*-

from bokeh.layouts import column
from bokeh.models import Div

class series_plot_cls(object) :
    do_gene_plot = True
    gene_plot = None
    do_selection_range = None
    do_stats_plots = False
    stats_plots = None
    stats_plot_names = [ 'samples_in_series', 'series_count', 'mean_series_length']
    match_test_allele_mask = None
    do_series_data_table = True
    series_data_table_html = None
    series_sample_match_data = None
    match_test_allele_mask = None
    series_plot_data = None
    plot_series_height_from_series_allele_count = True
    series_plot = None
    layout_items = None
    
    def assemble_layout(self) :
        self.layout_items = []
        if self.do_gene_plot or self.do_stats_plots :
            plots = []
            if self.gene_plot is not None :
                plots.append(self.gene_plot)
            if self.stats_plots is not None :
                plots.extend(self.stats_plots)
            plots.append(self.series_plot)
            self.series_column = column(plots)
        else :
            self.series_column = self.series_plot
        self.layout_items.append([self.series_column])
        if self.series_data_table_html is not None :
            self.layout_items.append([Div(text=self.series_data_table_html)])

class all_series_plot_cls(series_plot_cls) :
    do_stats_plots = True
    do_series_data_table = False


class select_series_plot_cls(series_plot_cls) :
    def __init__(self, select_series_index_and_id=None, select_allele_mask_and_id=None) :
        self.select_series_index_and_id = select_series_index_and_id
        self.select_allele_mask_and_id = select_allele_mask_and_id
        self.do_selection_range = True
        
class hierarchy_series_plot_cls(series_plot_cls) :
    def __init__(self, hierarchy_series_index_and_id, min_match=0.9):
        self.hierarchy_series_index_and_id = hierarchy_series_index_and_id
        self.min_match = min_match
        
class yes_no_series_plot_cls(series_plot_cls) :
    def __init__(self, yes_series_data_indexes_and_ids, no_series_data_indexes_and_ids=None, min_match=0.9) :
        self.yes_series_data_indexes_and_ids = yes_series_data_indexes_and_ids
        self.no_series_data_indexes_and_ids = no_series_data_indexes_and_ids
        self.min_match = min_match
        
class allele_mask_yes_no_series_plot_cls(yes_no_series_plot_cls) :
    def __init__(self, allele_mask_and_id, yes_series_data_indexes_and_ids=None,
                 no_series_data_indexes_and_ids=None, min_match=0.9) :
        yes_no_series_plot_cls.__init__(self, yes_series_data_indexes_and_ids, no_series_data_indexes_and_ids, min_match)
        self.yes_allele_mask_and_id = allele_mask_and_id
        