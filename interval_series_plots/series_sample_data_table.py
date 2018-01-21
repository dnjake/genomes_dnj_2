# -*- coding: utf-8 -*-
from __future__ import division

from ..html_display import array_table as html
from . import series_sample_data as ssd

num_tag = '<td style="text-align: right;">'
center_tag = '<td style="text-align: center;">'
int_fmt =  '{:d}'
big_fmt =  '{:,}'
float_fmt = '{:.2f}'
exp_fmt = '{:.1e}'

def region_stats_from_allele_masks(allele_masks) :
    series_stats = ssd.get_region_stats_from_allele_masks(allele_masks)
    r_c, o_p = series_stats
    cols = ((('afr',2), ((r_c['afr'], num_tag, int_fmt), (o_p['afr'], num_tag, float_fmt))),
           (('afx',2), ((r_c['afx'], num_tag, int_fmt), (o_p['afx'], num_tag, float_fmt))),
           (('amr',2), ((r_c['amr'], num_tag, int_fmt), (o_p['amr'], num_tag, float_fmt))),
           (('eas',2), ((r_c['eas'], num_tag, int_fmt), (o_p['eas'], num_tag, float_fmt))), 
           (('eur',2), ((r_c['eur'], num_tag, int_fmt), (o_p['eur'], num_tag, float_fmt))),
           (('sas',2), ((r_c['sas'], num_tag, int_fmt), (o_p['sas'], num_tag, float_fmt))),
           (('sax',2), ((r_c['sax'], num_tag, int_fmt), (o_p['sax'], num_tag, float_fmt))))
    return cols                       



class series_data_table_cls(object) :
    def __init__(self, po) :
        self.plot_data = po.plot_data
        self.match_test_allele_mask = po.match_test_allele_mask
        self.with_matches = True
        if self.plot_data['match_allele_count'] is None :
            self.with_matches = False

    def sort_by_allele_counts(self, allele_count) :
        arg_order = allele_count.argsort()
        arg_order = arg_order[::-1]
        self.sorted_plot_data = {}
        for key in self.plot_data.keys():
            data = self.plot_data[key]
            if data is None :
                self.sorted_plot_data[key] = None
            else :
                self.sorted_plot_data[key] = data[arg_order]


    def basic_series_data(self) :
        pd = self.sorted_plot_data
        cols = [('index', pd['data_index'], num_tag, int_fmt),
                ('first', pd['first_pos'], num_tag, big_fmt),
                ('last', pd['last_pos'], num_tag, big_fmt), 
                ('length', pd['series_length'], num_tag, big_fmt),
                ('snps', pd['snp_count'], num_tag, int_fmt)]
        if self.with_matches :
            match_counts_f = pd['match_allele_count'].astype('f4')
            allele_ratios = match_ratios = match_counts_f/pd['allele_count'].astype('f4')
            if self.match_test_allele_mask is None :
                alleles_and_matches = [
                    ('alleles', pd['allele_count'], num_tag, int_fmt), 
                    (('matches',2), ((pd['match_allele_count'], num_tag, int_fmt), (allele_ratios, num_tag, float_fmt)))
                ]
                cols.extend(alleles_and_matches)
            else :
                match_test_allele_count = float(self.match_test_allele_mask.sum())
                match_ratios = match_counts_f/match_test_allele_count
                alleles_and_matches = [
                    (('alleles',2), ((pd['allele_count'], num_tag, int_fmt), (allele_ratios, num_tag, float_fmt))),
                    (('matches',2), ((pd['match_allele_count'], num_tag, int_fmt), (match_ratios, num_tag, float_fmt)))
                ]
                cols.extend(alleles_and_matches)
        else :
             cols.append(('alleles', pd['allele_count'], num_tag, int_fmt))
             
        return cols
    
    def series_data_html(self) :
        self.sort_by_allele_counts(self.plot_data['allele_count'])
        data_cols  = self.basic_series_data()
        if self.with_matches :
            allele_masks = self.sorted_plot_data['match_allele_mask']
        else :
            allele_masks = self.sorted_plot_data['allele_mask']
        data_cols.extend(region_stats_from_allele_masks(allele_masks))
        html_table = html.html_table_cls(data_cols)
        return html_table.assemble_table()
    




















