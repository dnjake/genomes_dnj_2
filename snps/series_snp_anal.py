# -*- coding: utf-8 -*-
import numpy as np
from bokeh.models import Div
from . import series_snp_data as ssd
from . import series_snps_plot as ssp
from ..interval_series_plots import plot_style as style
from ..interval_series_plots import series_plot as sp
from ..interval_series_plots import series_anal_plots as sap

class series_snps_anal_cls(object) :
    series_snps_obj = None
    snps_per_sample = None
    selected_by_snp_count = None
    samples_per_snp = None
    yes_series_data_indexes_and_ids = None
    no_series_data_indexes_and_ids = None
    match_test_allele_mask = None
    min_match = 0.9
    plot_style = style
    table_style = ['<style type="text/css">']
    table_style.append('''td, th {
      border: 1px solid rgb(190,190,190);
      padding: 10px 20px;
    }''')
    table_style.append('</style>')
    table_style = '\n'.join(table_style)
    max_cols_per_row = 15
    def __init__(self, series_snps_obj, sample_select_mask=None) :
        self.series_snps_obj = series_snps_obj
        self.sample_select_mask = sample_select_mask
        self.series_id = self.series_snps_obj.series_id
        self.layout_items = []
    
    def table_header_row_html(self, header, data_array) :
        row = ['<tr><th scope=row style="text-align:right;padding-right:10px;">' + header + '</th>']
        for num in data_array :
            row.append('<th scope=col style="text-align:center;padding:0 10px;">' + str(num) + '</th>')
        row.append('</tr>')
        return ''.join(row)

    def table_row_html(self, header, data_array) :
        row = ['<tr><th scope=row style="text-align:right;padding-right:10px;">' + header + '</th>']
        for num in data_array :
            row.append('<td style="text-align:center;padding:0 10px">' + str(num) + '</td>')
        row.append('</tr>')
        return ''.join(row)

    def two_array_html(self, top_header, top_array, bottom_header, bottom_array) :
        html_table = [self.plot_style.table_tag]
        html_table.append(self.table_header_row_html(top_header, top_array))
        html_table.append(self.table_row_html(bottom_header, bottom_array))
        html_table.append('</table>')
        return '\n'.join(html_table)        
    
    def multi_table_html(self, top_header, top_array, bottom_header, bottom_array) :
        tables_html = []
        next_top_array = top_array
        next_bottom_array = bottom_array
        while next_top_array is not None :
            if next_top_array.size > self.max_cols_per_row :
                working_top = next_top_array[:self.max_cols_per_row]
                working_bottom = next_bottom_array[:self.max_cols_per_row]
                next_top_array = next_top_array[self.max_cols_per_row:]
                next_bottom_array = next_bottom_array[self.max_cols_per_row:]
            else :
                working_top = next_top_array
                working_bottom = next_bottom_array
                next_top_array = None
                next_bottom_array = None
            tables_html.append(self.two_array_html(top_header, working_top, bottom_header, working_bottom))
            tables_html.append('<p>')
        return '\n'.join(tables_html)

    def snps_per_sample_html(self) :
        if self.snps_per_sample is None :
            return
        top_header = 'snps'
        bottom_header = 'samples'
        top_array = self.snps_per_sample['count']
        bottom_array = self.snps_per_sample['snps']
        return self.multi_table_html(top_header, top_array, bottom_header, bottom_array)


    def samples_per_snp_selected_by_snp_count_html(self) :
        if self.selected_by_snp_count is None :
            return
        snp_count, samples_per_snp, snp_count_sample_mask = self.selected_by_snp_count
        snp_number = np.arange(samples_per_snp.size, dtype='i4')
        top_header = 'snp'
        bottom_header = 'samples'
        html = self.multi_table_html(top_header, snp_number, bottom_header, samples_per_snp)
        html = '<h4 style="min-width:500px">sample count by series snp</h4>\n' + html
        return html
        
    def series_alleles_per_snp_html(self) :
        alleles_per_snp = self.series_snps_obj.series_alleles_per_snp()
        snp_number = np.arange(alleles_per_snp.size, dtype='i4')
        top_header = 'snp'
        bottom_header = 'samples'
        html = self.multi_table_html(top_header, snp_number, bottom_header, alleles_per_snp)
        hd_html = ['<h4 style="min-width:500px">sample count by series ']
        hd_html.append(self.series_id)
        hd_html.append(' snp</h4>\n')
        html = ''.join(hd_html) + html
        return html
        
    def layout_append_series_aps_html(self):
        html = self.series_alleles_per_snp_html()
        self.layout_items.append([Div(text=html)])
    
    def match_test_samples_per_snp_html(self) :
        samples_per_snp = self.series_snps_obj.alleles_per_snp(self.match_test_allele_mask)
        snp_number = np.arange(samples_per_snp.size, dtype='i4')
        top_header = 'snp'
        bottom_header = 'samples'
        html = self.multi_table_html(top_header, snp_number, bottom_header, samples_per_snp)
        html = '<h4 style="min-width:500px">sample count by series snp</h4>\n' + html
        return html
        

class do_series_snps_anal_cls(object) :
    def __init__(self, series_anal_obj) :
        self.series_anal_obj = series_anal_obj    
        self.chrom = series_anal_obj.chrom
        self.anal_pos_range = self.series_anal_obj.anal_first_pos, self.series_anal_obj.anal_last_pos
        self.anal_selection_range = self.series_anal_obj.selection_range
        self.anal_series_plots = sap.series_anal_plt_cls(self.series_anal_obj)
        
    def create_anal_obj(self, series_data_index, sample_select_mask=None) :
        ai = self.series_anal_obj.array_index_from_data_index(series_data_index)
        sd = self.series_anal_obj.series_data[ai]
        am = self.series_anal_obj.allele_masks[ai]
        snps_obj = ssd.series_snps_cls.snps_for_series(self.chrom, sd, am)
        return series_snps_anal_cls(snps_obj, sample_select_mask)
    
    def do_snps_per_sample(self, ssao) :
        ssao.snps_per_sample = ssao.series_snps_obj.unique_snps_per_allele(ssao.sample_select_mask)
        return ssao

    def do_snps_per_sample_plot(self, ssao) :
        plt_obj = ssp.snp_match_count_plot_cls(ssao.snps_per_sample)
        snps_per_sample_plot = plt_obj.do_plot()
        ssao.layout_items.append([snps_per_sample_plot])

    def do_sps_header(self, ssao) :
        html = ['<h4 style="min-width:500px;">series ']
        html.append(ssao.series_id)
        html.append(' snp counts for ')
        html.append(ssao.sample_id)
        html.append(' samples</h4>\n')
        html = ''.join(html)
        ssao.layout_items.append([Div(text=html)])

    def do_sps_with_html(self, ssao) :
        self.do_snps_per_sample(ssao)
        self.do_sps_header(ssao)
        html = ssao.snps_per_sample_html()
        ssao.layout_items.append([Div(text=html)])
        
    def do_sps_with_plot_and_html(self, ssao) :
        self.do_snps_per_sample(ssao)
        self.do_sps_header(ssao)
        self.do_snps_per_sample_plot(ssao)
        html = ssao.snps_per_sample_html()
        ssao.layout_items.append([Div(text=html)])
        
    def do_select_samples_by_snp_count(self, ssao, snp_count) :
        aps, am = ssao.series_snps_obj.snps_from_aps_value(snp_count, ssao.sample_select_mask)
        ssao.selected_by_snp_count = snp_count, aps, am
        return ssao
    
    def do_select_by_snp_count_with_html(self, ssao, snp_count) :
        self.do_select_samples_by_snp_count(ssao, snp_count)
        html = []
        html.append('<h4 style="min-width:500px">')
        html.append(ssao.sample_id)
        html.append(' samples selected for ')
        html.append(str(snp_count))
        html.append(' ')
        html.append(ssao.series_id)
        html.append(' snps</h4>')
        html = ''.join(html)
        html = html + ssao.samples_per_snp_selected_by_snp_count_html()
        ssao.layout_items.append([Div(text=html)])
        
    def do_select_by_snp_count_with_html_and_plot(self, ssao, snp_count) :
        self.do_select_by_snp_count_with_html(ssao, snp_count)
        self.do_selected_allele_mask_series_plot(ssao)


    def do_samples_per_snp(self, ssao) :
        ssao.samples_per_snp = ssao.series_snps_obj.alleles_per_snp(ssao.sample_select_mask)
        return ssao
        
    def do_snp_count_series_plot(self, ssao):
        snp_count, snp_count_alleles_per_snp, snp_count_allele_mask = ssao.selected_by_snp_count
        id = 'am_' + ssao.sample_id + '_' + str(snp_count)
        am_and_id = snp_count_allele_mask, id        
        am_series_plot = sp.allele_mask_yes_no_series_plot_cls(am_and_id, ssao.yes_series_data_indexes_and_ids,
                                                               ssao.no_series_data_indexes_and_ids, ssao.min_match)
        self.anal_series_plots.do_allele_mask_yes_no_series_plot(am_series_plot)
        html = '<h4 style="min-width:500px;">snp count selected samples yes no plot</h4>\n'
        ssao.layout_items.append([Div(text=html)])
        ssao.layout_items.extend(am_series_plot.layout_items)
        ssao.match_test_allele_mask = am_series_plot.match_test_allele_mask
        
        
    def do_snp_anal_series_plot(self, ssao) :
        if ssao.sample_select_mask is None :
            series_plot = sp.yes_no_series_plot_cls(ssao.yes_series_data_indexes_and_ids,
                                                       ssao.no_series_data_indexes_and_ids, ssao.min_match)
            self.anal_series_plots.do_yes_no_series_plot(series_plot)
        else :
            id = 'am_' + ssao.sample_id
            am_and_id = ssao.sample_select_mask, id        
            series_plot = sp.allele_mask_yes_no_series_plot_cls(am_and_id, ssao.yes_series_data_indexes_and_ids,
                                                                   ssao.no_series_data_indexes_and_ids, ssao.min_match)
            self.anal_series_plots.do_allele_mask_yes_no_series_plot(series_plot)
        html = '<h4 style="min-width:500px;">snp anal ' + ssao.sample_id + ' yes no plot</h4>\n'
        ssao.layout_items.append([Div(text=html)])
        ssao.layout_items.extend(series_plot.layout_items)
        ssao.match_test_allele_mask = series_plot.match_test_allele_mask
    
    def do_match_test_samples_per_snp_html(self, ssao) :
        html = '<h4 style="min-width:500px;">plot match samples per series snp</h4>\n'
        html = html + ssao.match_test_samples_per_snp_html()
        ssao.layout_items.append(Div(text=html))
        
    def do_snp_count_plot_with_sps_html(self, ssao, snp_count=None) :
        if snp_count is not None :
            self.do_select_samples_by_snp_count(snp_count)
        self.do_snp_count_series_plot(ssao)
        self.do_match_test_samples_per_snp_html(ssao)
        
    def do_snp_anal_plot_with_sps_html(self, ssao) :
        self.do_snp_anal_series_plot(ssao)
        self.do_match_test_samples_per_snp_html(ssao)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        