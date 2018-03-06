# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
from itertools import izip

from . import series_sample_data_table as sdt
from ..html_display import array_table as html

class series_hierarchy_data_cls(object) :
    match_test_allele_mask = None
    series_sample_match_data = None
    series_plot_data = None
    root_series_data_indexes = None
    def __init__(self) :
        self.root_series_array_indexes = []
    

class series_hierarchy_anal_cls(object) :
    single_parent_min = 0.9
    child_parent_dtype = np.dtype([('child_data_index', 'u4'), ('child_snp_count', 'u2'), ('child_sample_count', 'u2'),
                                   ('parent_data_index', 'u4'),('parent_snp_count', 'u2'), ('parent_sample_count', 'u2'),
                                   ('match_allele_mask', '?', 5008)])
    def __init__(self, series_anal_obj) :
        #self.interval_obj = interval_obj
        self.series_anal_obj = series_anal_obj
        self.chrom = self.series_anal_obj.chrom
        #self.anal_pos_range = interval_obj.anal_pos_range
        self.series_data = self.series_anal_obj.series_data
        self.sample_data = self.series_anal_obj.sample_data
        self.allele_masks = self.series_anal_obj.allele_masks
        self.data_indexes = self.series_data['data_index']
        self.alleles_per_series = self.series_anal_obj.alleles_per_series
        self.series_order = self.alleles_per_series.argsort()
        self.single_child_parent_array_indexes = []
        self.root_series_array_indexes = []
                
    def look_for_parents(self, child_array_index, candidate_parent_array_indexes) :
        child_allele_mask = self.allele_masks[child_array_index]
        child_sample_count = float(child_allele_mask.sum())
        for ind in candidate_parent_array_indexes :
            candidate_allele_mask = self.allele_masks[ind]
            matches = np.logical_and(child_allele_mask, candidate_allele_mask)
            match_count = float(matches.sum())
            match_ratio = match_count/child_sample_count
            if match_ratio >= self.single_parent_min :
                self.single_child_parent_array_indexes.append((child_array_index, ind))
                return
        self.root_series_array_indexes.append(child_array_index)
        
    def package_root_data(self, hdo) :
        if len(self.root_series_array_indexes) == 0 :
            hdo.root_series_sample_data = None
            hdo.root_series_data_indexes = None
            return
        root_series_array_indexes = np.array(self.root_series_array_indexes)
        hdo.root_series_array_indexes = root_series_array_indexes[::-1]
        hdo.root_series_data = self.series_data[hdo.root_series_array_indexes]
        hdo.root_sample_data = self.sample_data[hdo.root_series_array_indexes]
        hdo.root_series_sample_data = {'series_data': hdo.root_series_data,
                                        'sample_data': hdo.root_sample_data,
                                        'match_data': None}
        hdo.root_series_data_indexes = hdo.root_series_data['data_index'].copy()
        hdo.root_series_snp_counts = self.series_data['item_count'][hdo.root_series_array_indexes]
        hdo.root_series_sample_counts = self.alleles_per_series[hdo.root_series_array_indexes]
        hdo.root_series_data_indexes_and_ids = []
        for data_index, snp_count, sample_count in izip(hdo.root_series_data_indexes, 
                                                        hdo.root_series_snp_counts, hdo.root_series_sample_counts) :
            hdo.root_series_data_indexes_and_ids.append((data_index, str(snp_count) + '_' + str(sample_count)))
            
            
    def package_child_parent_data(self, hdo) :
        if len(self.single_child_parent_array_indexes) == 0 :
            hdo.child_parent_data = None
            return
        child_array_indexes, parent_array_indexes = zip(*self.single_child_parent_array_indexes)
        child_array_indexes = np.array(child_array_indexes)
        parent_array_indexes = np.array(parent_array_indexes)
        cpd = np.zeros(child_array_indexes.size, self.child_parent_dtype)
        csd = self.series_data[child_array_indexes]
        psd = self.series_data[parent_array_indexes]
        cpd['child_data_index'] = csd['data_index']
        cpd['child_snp_count'] = csd['item_count']
        cpd['child_sample_count'] = csd['p90_allele_count']
        cpd['parent_data_index'] = psd['data_index']
        cpd['parent_snp_count'] = psd['item_count']
        cpd['parent_sample_count'] = psd['p90_allele_count']
        cpd['match_allele_mask'] = np.logical_and(self.allele_masks[child_array_indexes], self.allele_masks[parent_array_indexes])
        hdo.child_parent_data = cpd[::-1]
        
    def get_root_data_html(self, hdo) :
        if hdo.root_series_sample_data is None :
            self.root_data_table_html = None
            return
    
        hdo.plot_data = self.series_anal_obj.plot_data_from_series_sample_match_data(hdo.root_series_sample_data)
        #self.match_test_allele_mask = None
        rdt_obj = sdt.series_data_table_cls(hdo)
        hdo.root_series_data_html = rdt_obj.series_data_html()
        
    def get_child_parent_html(self, hdo) :
        if hdo.child_parent_data is None :
            hdo.child_parent_data_table_html = None
            return
        cpd = hdo.child_parent_data
        cdi = cpd['child_data_index']
        csnpc = cpd['child_snp_count']
        csampc = cpd['child_sample_count']
        pdi = cpd['parent_data_index']
        psnpc = cpd['parent_snp_count']
        psampc = cpd['parent_sample_count']
        cols = [('child index', cdi, sdt.num_tag, sdt.int_fmt),
                ('child snps', csnpc, sdt.num_tag, sdt.int_fmt),
                ('child alleles', csampc, sdt.num_tag, sdt.int_fmt),
                ('parent index', pdi, sdt.num_tag, sdt.int_fmt),
                ('parent snps', psnpc, sdt.num_tag, sdt.int_fmt),
                ('parent alleles', psampc, sdt.num_tag, sdt.int_fmt) ]
        match_count = hdo.child_parent_data['match_allele_mask'].sum(axis=1)
        match_count_f = match_count.astype('f4')
        samp_count_f = csampc.astype('f4')
        match_ratios = match_count_f/samp_count_f
        matches = (('matches',2), ((match_count, sdt.num_tag, sdt.int_fmt), (match_ratios, sdt.num_tag, sdt.float_fmt)))
        cols.append(matches)
        cols.extend(sdt.region_stats_from_allele_masks(hdo.child_parent_data['match_allele_mask']))
        html_table = html.html_table_cls(cols)
        hdo.child_parent_series_data_html = html_table.assemble_table()
    
    def find_hierarchy_roots(self, hdo) :    
        for ind in range(self.series_order.size) :
            child_array_index = self.series_order[ind]
            start_rest = ind + 1
            candidate_parent_array_indexes = self.series_order[start_rest:]
            self.look_for_parents(child_array_index, candidate_parent_array_indexes)
        self.package_root_data(hdo)
        self.package_child_parent_data(hdo)
        
    def do_hierarchy(self) :
        self.hierarchy_data_obj = series_hierarchy_data_cls()
        self.find_hierarchy_roots(self.hierarchy_data_obj)
        self.get_root_data_html(self.hierarchy_data_obj)
        self.get_child_parent_html(self.hierarchy_data_obj)
        return self.hierarchy_data_obj
        
