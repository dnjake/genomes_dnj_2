# -*- coding: utf-8 -*-
import numpy as np
from itertools import izip
from bokeh.plotting import output_file, show
from bokeh.layouts import layout, column
from bokeh.models import Div
from . import series_hierarchy_anal as sha
from . import series_anal_plots as sap
from .series_plot import all_series_plot_cls as aspc, select_series_plot_cls as sspc, \
                                                         hierarchy_series_plot_cls as hspc
from . import plot_style as style
class standard_page_cls(object) :
    plot_style = style
    tab_title = None

    def __init__(self) :
        self.page_layout_items = []

    def show_page(self) :
        bottom_html = self.plot_style.bottom_margin_div()
        self.page_layout_items.append([Div(text=bottom_html)])
        if self.tab_title is None :
            output_file(self.html_file_name)
        else :
            output_file(self.html_file_name, title=self.tab_title)
        show(layout(self.page_layout_items))
    
    
class all_series_page_cls(standard_page_cls) :
    tab_title = 'all_interval_series'
    def __init__(self, do_series_plots_obj, hierarchy_data_obj) :
        standard_page_cls.__init__(self)
        self.do_series_plots_obj = do_series_plots_obj
        self.hierarchy_data_obj = hierarchy_data_obj
        interval_id = self.do_series_plots_obj.interval_id
        self.all_series_id = 'all_series_' + interval_id
        self.html_file_name = self.all_series_id + '.html'
        self.tab_title = self.all_series_id
     
    def do_plot(self) :
        po = aspc()
        po.do_series_data_table = True
        po.do_selection_range = True
        self.do_series_plots_obj.do_all_series_plot(po)
        self.page_layout_items.extend(po.layout_items)

    def do_hierarchy_html(self) :
        root_html = self.hierarchy_data_obj.root_series_data_html
        if root_html is not None :
            hd = '<h4>Series Hierarchy Root Series</h4>\n'
            root_html = hd + root_html
            self.page_layout_items.append([Div(text=root_html)])
        child_parent_html = self.hierarchy_data_obj.child_parent_series_data_html
        if child_parent_html is not None :
            hd = '<h4>Series Child To Parent Relations</h4>\n'
            child_parent_html = hd + child_parent_html
            self.page_layout_items.append(Div(text=child_parent_html))
        
    def do_page(self) :
        self.do_plot()
        self.do_hierarchy_html()
        self.show_page()

class selected_series_page_cls(standard_page_cls) :
    tab_title = 'selected_interval_series'
    def __init__(self, do_series_plots_obj) :
        standard_page_cls.__init__(self)
        self.do_series_plots_obj = do_series_plots_obj
        self.interval_id = self.do_series_plots_obj.interval_id
        self.selected_id = 'selected_series_' + self.interval_id
        self.html_file_name = self.selected_id + '.html'
        self.tab_title = self.selected_id

    def do_most_common_with_stats(self) :
        sspo = sspc()
        sspo.do_stats_plots = True
        self.do_series_plots_obj.do_most_common_series_plot(sspo)
        self.page_layout_items.extend(sspo.layout_items)
        
    def do_basis(self) :
        #plot_title = 'basis_series_' + self.interval_id
        bspo = self.do_series_plots_obj.do_basis_series_plot(sspc())
        self.page_layout_items.extend(bspo.layout_items)
        
    def do_longest(self) :        
        #plot_title = 'longest_series_' + self.interval_id
        lspo = self.do_series_plots_obj.do_longest_series_plot(sspc())
        self.page_layout_items.extend(lspo.layout_items)
        
    def do_page(self) :
        self.do_most_common_with_stats()
        self.do_basis()
        self.do_longest()
        self.show_page()

class series_hierarchy_page_cls(standard_page_cls) :
    high_sample_data_match_min = 0.9
    low_sample_data_match_min = 0.1
    def __init__(self, do_series_plots_obj, series_data_index_and_id) :
        standard_page_cls.__init__(self)
        self.do_series_plots_obj = do_series_plots_obj
        self.series_data_index_and_id = series_data_index_and_id
        self.series_data_index, self.series_id = self.series_data_index_and_id
        self.series_id = self.series_id + '_' + str(self.series_data_index)
        self.html_file_name = 'sh_' + self.series_id + '.html'
        self.tab_title = self.series_id
    
    def do_hierarchy_plot(self) :    
        #plot_title = 'high_match_' + self.series_id
        high_hspo = hspc(self.series_data_index_and_id, min_match=self.high_sample_data_match_min)
        high_hspo.do_selection_range = True
        self.do_series_plots_obj.do_hierarchy_series_plot(high_hspo)
        self.page_layout_items.extend(high_hspo.layout_items)
        low_hspo = hspc(self.series_data_index_and_id, min_match=self.low_sample_data_match_min)
        low_hspo.do_selection_range = True
        self.do_series_plots_obj.do_hierarchy_series_plot(low_hspo)
        self.page_layout_items.extend(low_hspo.layout_items)
    
    def do_basis(self) :
        bspo = sspc(self.series_data_index_and_id)
        self.do_series_plots_obj.do_basis_series_plot(bspo)
        self.page_layout_items.extend(bspo.layout_items)
        
    def do_page(self) :
        self.do_hierarchy_plot()
        self.do_basis()
        self.show_page()

class interval_standard_plots_cls(object) :
    max_root_series = 20
    min_longest_series_matches = 128
    #series_stats_names = [ 'samples_in_series', 'series_count', 'mean_series_length']

    def __init__(self, do_series_plots_obj) :
        self.do_series_plots_obj = do_series_plots_obj
        self.series_anal_obj = self.do_series_plots_obj.series_anal_obj

    @classmethod
    def all_read_series_from_chrom_data(cls, chrom, first_pos, last_pos) :
        sapc = sap.series_anal_plt_cls
        sapo = sapc.all_read_series_from_chrom_data(chrom, first_pos, last_pos)
        return cls(sapo)

    @classmethod 
    def selected_series_from_chrom_data(cls, chrom, first_pos, last_pos) :
        sapc = sap.series_anal_plt_cls
        sapo = sapc.selected_series_from_chrom_data(chrom, first_pos, last_pos)
        return cls(sapo)
    
    @classmethod
    def selected_series_from_input_data(cls, in_data_obj, first_pos, last_pos) :
        sapc = sap.series_anal_plt_cls
        sapo = sapc.selected_series_from_input_data(in_data_obj, first_pos, last_pos)
        return cls(sapo)

    def find_hierarchy_roots(self) :
        hierarchy_obj = sha.series_hierarchy_anal_cls(self.series_anal_obj)
        self.hierarchy_data_obj = hierarchy_obj.do_hierarchy()

    def find_longest_series(self) :
        self.do_series_plots_obj.selection_obj.select_by_length()
        longest_series_sample_data = self.do_series_plots_obj.selection_obj.select_by_length()
        lsm = longest_series_sample_data['match_data']
        match_counts = lsm['match_allele_mask'].sum(axis=1)
        m = match_counts >= self.min_longest_series_matches
        sd = longest_series_sample_data['series_data']
        data_indexes = sd['data_index'][m]
        snp_counts = sd['item_count'][m]
        sample_counts = sd['p90_allele_count'][m]
        self.longest_data_indexes_and_ids = []
        for data_index, snp_count, sample_count in izip(data_indexes, snp_counts, sample_counts) :
            self.longest_data_indexes_and_ids.append((data_index, str(snp_count) + '_' + str(sample_count)))
        

    def do_all_series(self) :
        page = all_series_page_cls(self.do_series_plots_obj, self.hierarchy_data_obj)
        page.do_page()

    def do_hierarchy_series_plots(self, series_data_indexes_and_ids) :
        for did in series_data_indexes_and_ids :
            po = series_hierarchy_page_cls(self.do_series_plots_obj, did)
            po.do_page()
        
    def do_root_series(self) :
        plot_series_dids = self.hierarchy_data_obj.root_series_data_indexes_and_ids
        if len(plot_series_dids) > self.max_root_series :
            plot_series_dids = plot_series_dids[:self.max_root_series]
        plot_data_indexes = set()
        for data_index, id in plot_series_dids :
            plot_data_indexes.add(data_index)
        for did in self.longest_data_indexes_and_ids :
            di, id = did
            if di not in plot_data_indexes :
                plot_series_dids.append(did)
                plot_data_indexes.add(di)
        self.do_hierarchy_series_plots(plot_series_dids)

    def do_selected_series(self) :
        sso = selected_series_page_cls(self.do_series_plots_obj)
        sso.do_page()

    
    def do_standard(self) :    
        self.find_hierarchy_roots()
        self.find_longest_series()
        self.do_all_series()
        self.do_selected_series()
        self.do_root_series()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    