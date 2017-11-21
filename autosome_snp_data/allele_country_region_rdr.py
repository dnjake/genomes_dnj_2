# -*- coding: utf-8 -*-

import os
import numpy as np
import tables as tb
from ..html_display import array_table as html

s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

in_file = 'allele_country_and_region_plus_maps.h5'
in_path = os.path.join(mod_dir, in_file)

class country_region_alleles_cls(object) :
    simple_region_codes = np.array(['afr', 'amr', 'eas', 'eur', 'sas'], dtype='S3')
    region_count_dtype = np.dtype([('afr', 'u2'), ('afx', 'u2'), ('amr', 'u2'), ('eas', 'u2'), 
                                   ('eur', 'u2'), ('sas', 'u2'), ('sax', 'u2')])                                   
    obs_to_pred_dtype = np.dtype([('afr', 'f4'), ('afx', 'f4'), ('amr', 'f4'), ('eas', 'f4'), 
                                   ('eur', 'f4'), ('sas', 'f4'), ('sax', 'f4')])                                   
    analysis_result_descr = [('code', 'S3'), ('allele_count', '<f4'), ('freq', '<f4'), ('obs_to_pred', '<f4')]
    analysis_result_dtype = np.dtype(analysis_result_descr)
    total_alleles = 5008
    float_total_alleles = float(total_alleles)    
    country_html_cols = 5
    
    def __init__(self) : 
        h5 = tb.open_file(in_path, 'r')
        self.countries = h5.root.country_alleles[:]
        self.regions = h5.root.region_alleles[:]
        h5.close()
        self.alleles_per_country = self.countries['allele_mask'].sum(axis=1)
        self.alleles_per_country = self.alleles_per_country.astype('f4')
        self.freq_alleles_per_country = self.alleles_per_country/self.float_total_alleles
        self.alleles_per_region = self.regions['allele_mask'].sum(axis=1)
        self.alleles_per_region = self.alleles_per_region.astype('f4')
        self.simple_alleles_per_region = np.zeros(self.simple_region_codes.size, dtype='f4')
        self.simple_alleles_per_region[[1, 2, 3]] = self.alleles_per_region[[2, 3, 4]]
        self.simple_alleles_per_region[0] = self.alleles_per_region[0] + self.alleles_per_region[1]
        self.simple_alleles_per_region[4] = self.alleles_per_region[5] + self.alleles_per_region[6]
        self.freq_alleles_per_region = self.alleles_per_region/self.float_total_alleles
        self.freq_simple_alleles_per_region = self.simple_alleles_per_region/self.float_total_alleles
        self.country_codes = self.countries['code']
        self.region_codes = self.regions['code']
        self.calc_country_cols()
        
    def calc_country_cols(self) :
        self.country_html_bounds = []
        remaining_countries = self.countries.size
        remaining_cols = self.country_html_cols
        start = 0
        while start < self.countries.size :
            col_countries = remaining_countries/remaining_cols
            if remaining_countries%remaining_cols > 0 :
                col_countries += 1
            bound = start + col_countries
            self.country_html_bounds.append((start, bound))
            remaining_countries -= col_countries
            remaining_cols -= 1
            start = bound

    def analyze_countries_from_allele_mask(self, allele_mask) :
        sample_country_alleles = np.logical_and(self.countries['allele_mask'], allele_mask)
        results = np.zeros(self.countries.size, self.analysis_result_dtype)
        results['code'] = self.country_codes
        results['allele_count'] = sample_country_alleles.sum(axis=1)
        results['freq'] = results['allele_count']/self.alleles_per_country
        alleles_in_mask = allele_mask.sum()
        pred = float(alleles_in_mask)*self.freq_alleles_per_country
        results['obs_to_pred'] = results['allele_count']/pred
        return results

    def analyze_countries(self, allele_indexes) :
        allele_mask = np.zeros(self.total_alleles, '?')
        allele_mask[allele_indexes] = True
        return self.analyze_countries_from_allele_mask(allele_mask)
    
    def analyze_regions(self, allele_indexes) :
        sample_region_alleles = self.regions['allele_mask'][:, allele_indexes]
        results = np.zeros(self.regions.size, self.analysis_result_dtype)
        results['code'] = self.region_codes
        results['allele_count'] = sample_region_alleles.sum(axis=1)
        results['freq'] = results['allele_count']/self.alleles_per_region
        pred = float(allele_indexes.size)*self.freq_alleles_per_region
        m = pred > 0
        results['obs_to_pred'][m] = results['allele_count'][m]/pred[m]
        return results

    def simple_region_stats(self, allele_indexes) :
        trs = self.analyze_regions(allele_indexes)
        results = np.zeros(self.simple_region_codes.size, self.analysis_result_dtype)
        results['code'] = self.simple_region_codes
        trs_counts = trs['allele_count']
        results_counts = results['allele_count']
        results_counts[[1, 2, 3]] = trs_counts[[2, 3, 4]]
        results_counts[0] = trs_counts[0] + trs_counts[1]
        results_counts[4] = trs_counts[5] + trs_counts[6]
        results['freq'] = results['allele_count']/self.simple_alleles_per_region
        pred = float(allele_indexes.size)*self.freq_simple_alleles_per_region
        m = pred > 0
        results['obs_to_pred'][m] = results['allele_count'][m]/pred[m]
        return results

    def analyze_regions_and_countries(self, allele_indexes) :
        country_results = self.analyze_countries(allele_indexes)
        region_results = self.analyze_regions(allele_indexes)
        return region_results, country_results

    def region_allele_indexes(self, allele_indexes, region_code) :
        region_code = region_code.upper()
        region_index = self.region_codes.searchsorted(region_code)
        region_allele_indexes = np.where(self.regions['allele_mask'][region_index,:])[0]
        return np.intersect1d(allele_indexes, region_allele_indexes)
        
    def region_allele_mask(self, region_code) :
        region_code = region_code.upper()
        region_index = self.region_codes.searchsorted(region_code)
        return self.regions['allele_mask'][region_index,:].copy()

    def region_counts_from_allele_masks(self, allele_masks) :
        data_rows = allele_masks.shape[0]
        out_data = np.zeros(data_rows, self.region_count_dtype)
        for code in self.region_codes :
            code_allele_mask = self.region_allele_mask(code)
            and_allele_masks = np.logical_and(allele_masks, code_allele_mask)
            out_data[code.lower()] = and_allele_masks.sum(axis=1)
        return out_data

    def region_obs_to_pred_from_region_counts(self, region_counts, maybe_allele_mask=None) :
        if maybe_allele_mask is None :
            maybe_region_counts = self.alleles_per_region
        else :
            maybe_region_counts = self.maybe_region_counts_from_maybe_allele_mask(maybe_allele_mask)
            maybe_region_counts = maybe_region_counts.astype('f4')
        total_maybe_counts = maybe_region_counts.sum()
        obs_to_pred = []
        for ind_region_counts in region_counts :
            ind_region_counts = np.array(list(ind_region_counts), dtype='f4')
            ind_total_obs_counts = ind_region_counts.sum()
            pred_counts = (maybe_region_counts/total_maybe_counts)*ind_total_obs_counts
            m = pred_counts > 0.0
            ind_region_counts[m] = ind_region_counts[m]/pred_counts[m]
            obs_to_pred.append(tuple(ind_region_counts))
        obs_to_pred = np.array(obs_to_pred, dtype=self.obs_to_pred_dtype)
        return obs_to_pred
        
    def maybe_regions_from_maybe_allele_mask(self, maybe_allele_mask) :
        regions = self.regions.copy()
        region_allele_masks = regions['allele_mask']
        for ind in range(region_allele_masks.shape[0]) :
            region_allele_masks[ind] = np.logical_and(region_allele_masks[ind], maybe_allele_mask)
        return regions
        
    def maybe_region_counts_from_maybe_allele_mask(self, maybe_allele_mask) :
        maybe_regions = self.maybe_regions_from_maybe_allele_mask(maybe_allele_mask)
        maybe_alleles_per_region = maybe_regions['allele_mask'].sum(axis=1)
        return maybe_alleles_per_region
        
    def maybe_region_counts_from_allele_masks(self, allele_masks, maybe_allele_mask) :
        maybe_regions = self.maybe_regions_from_maybe_allele_mask(maybe_allele_mask)
        data_rows = allele_masks.shape[0]
        out_data = np.zeros(data_rows, self.region_count_dtype)
        for maybe_region in maybe_regions :
            code = maybe_region['code']
            code_allele_mask = maybe_region['allele_mask']
            and_allele_masks = np.logical_and(allele_masks, code_allele_mask)
            out_data[code.lower()] = and_allele_masks.sum(axis=1)
        return out_data

    def maybe_region_stats(self, allele_indexes, maybe_allele_mask) :
        maybe_alleles_total = float(maybe_allele_mask.sum())
        maybe_regions = self.maybe_regions_from_maybe_allele_mask(maybe_allele_mask)
        maybe_alleles_per_region = maybe_regions['allele_mask'].sum(axis=1)
        maybe_alleles_per_region = maybe_alleles_per_region.astype('f4')
        sample_region_alleles = maybe_regions['allele_mask'][:, allele_indexes]
        results = np.zeros(self.regions.size, self.analysis_result_dtype)
        results['code'] = self.region_codes
        allele_counts = sample_region_alleles.sum(axis=1)
        allele_counts = allele_counts.astype('f4')
        results['allele_count'] = allele_counts
        results['freq'] = results['allele_count']/maybe_alleles_per_region
        freq_maybe_alleles_per_region = maybe_alleles_per_region/maybe_alleles_total
        pred = float(allele_indexes.size)*freq_maybe_alleles_per_region
        m = pred > 0
        results['obs_to_pred'][m] = results['allele_count'][m]/pred[m]
        return results
        
    def maybe_simple_region_stats(self, allele_indexes, maybe_allele_mask) :
        maybe_alleles_total = float(maybe_allele_mask.sum())
        maybe_regions = self.maybe_regions_from_maybe_allele_mask(maybe_allele_mask)
        maybe_alleles_per_region = maybe_regions['allele_mask'].sum(axis=1)
        maybe_alleles_per_region = maybe_alleles_per_region.astype('f4')
        maybe_simple_alleles_per_region = np.zeros(self.simple_region_codes.size, dtype='f4')
        maybe_simple_alleles_per_region[[1, 2, 3]] = maybe_alleles_per_region[[2, 3, 4]]
        maybe_simple_alleles_per_region[0] = maybe_alleles_per_region[0] + maybe_alleles_per_region[1]
        maybe_simple_alleles_per_region[4] = maybe_alleles_per_region[5] + self.maybe_alleles_per_region[6]
        freq_maybe_simple_alleles_per_region = maybe_simple_alleles_per_region/maybe_alleles_total
        trs = self.maybe_region_stats(allele_indexes, maybe_allele_mask)
        results = np.zeros(self.simple_region_codes.size, self.analysis_result_dtype)
        results['code'] = self.simple_region_codes
        trs_counts = trs['allele_count']
        results_counts = results['allele_count']
        results_counts[[1, 2, 3]] = trs_counts[[2, 3, 4]]
        results_counts[0] = trs_counts[0] + trs_counts[1]
        results_counts[4] = trs_counts[5] + trs_counts[6]
        results['freq'] = results['allele_count']/maybe_simple_alleles_per_region
        pred = float(allele_indexes.size)*freq_maybe_simple_alleles_per_region
        m = pred > 0
        results['obs_to_pred'][m] = results['allele_count'][m]/pred[m]
        return results
        
    def data_html(self, data) :
        s_d = data
        allele_counts = s_d['allele_count'].astype('u2')
        num_tag = '<td style="text-align: right;">'
        str_tag = '<td style="text-align: center">'
        int_fmt =  '{:d}'
        float_fmt = '{:.2f}'
        column_info = (('pop', s_d['code'], str_tag),
                       ('alleles', allele_counts, num_tag, int_fmt),
                       ('freq', s_d['freq'], num_tag, float_fmt),
                       ('obs_to_pred', s_d['obs_to_pred'], num_tag, float_fmt)) 
        data_html_obj = html.html_table_cls(column_info)
        table_html = data_html_obj.assemble_table()
        return table_html

    def region_and_country_html(self, allele_indexes) :
        region_results, country_results = self.analyze_regions_and_countries(allele_indexes)
        out_str = ['<table style="width:800px;border-style:hidden;border-width:0;" >',
                    '<tr>',
                    '<td style="border-style:none;vertical-align:text-top;margin-left:100px;"><p>Allele Regions</p>', 
                    self.data_html(region_results), '</td>', 
                   '<td style="border-style:none;vertical-align:text-top;"><p>Allele Countries</p>',
                   self.data_html(country_results), '</td></tr></table>' ]
        return ''.join(out_str)
                   
    def country_html(self, country_data) :
        c_d = country_data
        allele_counts = c_d['allele_count'].astype('u2')
        num_tag = '<td style="text-align: right;">'
        str_tag = '<td style="text-align: center">'
        int_fmt =  '{:d}'
        float_fmt = '{:.2f}'
        column_info = (('pop', c_d['code'], str_tag),
                       ('cnt', allele_counts, num_tag, int_fmt),
                       ('otp', c_d['obs_to_pred'], num_tag, float_fmt)) 
        data_html_obj = html.html_table_cls(column_info)
        table_html = data_html_obj.assemble_table()
        return table_html
        
    def country_html_table(self, allele_mask) :
        country_data = self.analyze_countries_from_allele_mask(allele_mask)
        out_str = ['<table style="width:800px;border-style:hidden;border-width:0;" >''<tr>\n']
        for start, bound in self.country_html_bounds :
            out_str.append('<td style="border-style:none;vertical-align:text-top;">')
            out_str.append(self.country_html(country_data[start:bound]))
            out_str.append('</td>\n')
        out_str.append('</tr></table>\n')
        return ''.join(out_str)

    def country_html_table_from_allele_indexes(self, allele_indexes) :
        allele_mask = np.zeros(self.total_alleles, '?')
        allele_mask[allele_indexes] = True
        return self.country_html_table(allele_mask)