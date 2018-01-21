# -*- coding: utf-8 -*-
from __future__ import division

import numpy as np

class select_series_cls(object) :
    total_allele_count = 5008
    selection_data_dtype = np.dtype([('data_index', 'u4'), ('match_allele_mask', '?', total_allele_count)])

    def __init__(self, series_data_obj) :
        self.series_data_obj = series_data_obj
        #self.cra = self.interval_data.cra
        self.anal_first_pos = self.series_data_obj.anal_first_pos
        self.anal_last_pos = self.series_data_obj.anal_last_pos
        self.series_data = self.series_data_obj.series_data
        self.series_sample_data = self.series_data_obj.sample_data
        self.series_allele_masks = self.series_data_obj.allele_masks
        self.all_allele_mask = np.zeros(self.total_allele_count, '?')
        self.all_allele_mask[:] = True
        self.series_first_pos = self.series_data_obj.series_first_pos
        self.series_last_pos = self.series_data_obj.series_last_pos
        self.series_lengths = self.series_last_pos - self.series_first_pos + 1
        self.series_snp_counts = self.series_data['item_count']
        self.alleles_per_series = self.series_data_obj.alleles_per_series

    def to_out_data(self, selection_data) :
        selection_data.sort(order='data_index')
        out_data = self.series_data_obj.series_sample_data_from_data_indexes(selection_data['data_index'])
        out_data['match_data'] = selection_data
        return out_data

    def select_series_in_order(self, series_array_indexes, selection_allele_mask=None) :
        data_indexes = self.series_data['data_index']
        series_allele_masks = self.series_allele_masks
        #allele_counts = self.alleles_per_series
        not_selected_allele_mask = selection_allele_mask
        if not_selected_allele_mask is None :
            not_selected_allele_mask = self.all_allele_mask.copy()
        not_selected_count = not_selected_allele_mask.sum()
        selection_data = []
        for index in series_array_indexes :
            index_allele_mask = series_allele_masks[index]
            index_selected_mask = np.logical_and(not_selected_allele_mask, index_allele_mask)
            index_selected_count = index_selected_mask.sum()
            if index_selected_count == 0 :
                continue
            index_not_selected = np.logical_not(index_selected_mask)
            not_selected_allele_mask = np.logical_and(not_selected_allele_mask, index_not_selected)
            index_data = (data_indexes[index], index_selected_mask)
            selection_data.append(index_data)
            not_selected_count -= index_selected_count
            if not_selected_count == 0 :
                break                                             
        selection_data = np.array(selection_data, dtype=self.selection_data_dtype)
        return self.to_out_data(selection_data)
        
    def select_by_most_common(self) :
        sorted_indexes = self.alleles_per_series.argsort()
        sorted_indexes = sorted_indexes[::-1]
        return self.select_series_in_order(sorted_indexes)
        
    def select_by_snp_count(self) :
        snp_counts = self.series_snp_counts
        sorted_indexes = snp_counts.argsort()
        sorted_indexes = sorted_indexes[::-1]
        return self.select_series_in_order(sorted_indexes)
        
    def select_by_length(self) :
        first_pos = self.series_first_pos
        last_pos = self.series_last_pos
        lengths = last_pos - first_pos
        sorted_indexes = lengths.argsort()
        sorted_indexes = sorted_indexes[::-1]
        return self.select_series_in_order(sorted_indexes)
        
    def select_by_basis(self, selection_allele_mask=None) :
        sorted_indexes = self.alleles_per_series.argsort()
        return self.select_series_in_order(sorted_indexes, selection_allele_mask)

    def select_by_basis_from_plt_obj(self, po) :
        selection_allele_mask = None
        if po.select_series_data_index is not None :
            tdi = po.select_series_data_index
            selection_allele_mask = self.series_data_obj.allele_mask_from_data_index(tdi)
        elif po.select_allele_mask is not None :
            selection_allele_mask = po.select_allele_mask
        po.series_sample_match_data = self.select_by_basis(selection_allele_mask)

    def select_by_min_first_pos(self) :
        sorted_indexes = self.series_first_pos.argsort()
        return self.select_series_in_order(sorted_indexes)
        
    def select_by_max_last_pos(self) :
        sorted_indexes = self.series_last_pos.argsort()
        sorted_indexes = sorted_indexes[::-1]
        return self.select_series_in_order(sorted_indexes)
        
