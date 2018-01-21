# -*- coding: utf-8 -*-
from __future__ import division

import numpy as np

class series_data_anal_cls(object) :
    @staticmethod
    def plot_data_from_series_sample_match_data(series_sample_match_data) :
        plot_data = {}
        series_data = series_sample_match_data['series_data']
        sample_data = series_sample_match_data['sample_data']
        series_allele_masks = sample_data['allele_mask']
        match_data = series_sample_match_data['match_data']
        if match_data is None :
            match_allele_masks = None
            match_allele_counts = None
        else :
            match_allele_masks = match_data['match_allele_mask']
            match_allele_counts = match_allele_masks.sum(axis=1)
        for name in ['data_index', 'first_pos', 'last_pos'] :
            plot_data[name] = series_data[name]
        plot_data['series_length'] = series_data['last_pos'] - series_data['first_pos']
        plot_data['snp_count'] = series_data['item_count']
        plot_data['allele_count'] = series_allele_masks.sum(axis=1)
        plot_data['allele_mask'] = series_allele_masks
        plot_data['match_allele_count'] = match_allele_counts
        plot_data['match_allele_mask'] = match_allele_masks
        plot_data['series_data'] = series_data
        return plot_data

    sample_data_dtype = np.dtype([('data_index',  'u4'), ('allele_mask', '?', 5008)])
    match_data_dtype = np.dtype([('data_index',  'u4'), ('match_allele_mask', '?', 5008)])
    total_allele_count = 5008
    def __init__(self, chrom, series_sample_data, selection_range=None) :
        series_data = series_sample_data['series_data']
        sample_data = series_sample_data['sample_data']
        self.chrom = chrom
        self.series_data = series_data
        self.sample_data = sample_data
        self.allele_masks = self.sample_data['allele_mask']
        self.data_indexes = self.series_data['data_index']
        self.series_first_pos = self.series_data['first_pos']
        self.series_last_pos = self.series_data['last_pos']
        self.anal_first_pos = self.series_first_pos.min()
        self.anal_last_pos = self.series_last_pos.max()
        self.alleles_per_series = self.allele_masks.sum(axis=1)
        self.selection_range = selection_range
        
    def array_index_from_data_index(self, data_index) :
        return self.data_indexes.searchsorted(data_index)

    def series_data_from_data_indexes(self, data_indexes) :
        indexes = self.data_indexes.searchsorted(data_indexes)
        series_data = self.series_data[indexes]
        allele_masks = self.allele_masks[indexes]
        return series_data, allele_masks

    def simple_id_from_data_index(self, data_index) :
        array_index = self.data_indexes.searchsorted(data_index)
        sd = self.series_data[array_index]
        snp_count = sd['item_count']
        sample_count = sd['p90_allele_count']
        id = str(snp_count) + '_' + str(sample_count)
        return id
    
    def unique_id_from_data_index(self, data_index) :
        simple_id = self.simple_id_from_data_index(data_index)
        return simple_id + '_' + str(data_index)

    def series_sample_data_from_data_indexes(self, data_indexes) :
        indexes = self.data_indexes.searchsorted(data_indexes)
        series_data = self.series_data[indexes]
        sample_data = self.sample_data[indexes]
        out_data = {'series_data': series_data,
                    'sample_data': sample_data,
                    'match_data': None}
        return out_data

    def out_series_sample_data(self) :
        out_dict = {'series_data': self.series_data,
                    'sample_data': self.sample_data,
                    'match_data': None}
        return out_dict

    def out_data_from_data_mask(self, mask) :
        out_dict = {'series_data': self.series_data[mask],
                    'sample_data': self.sample_data[mask],
                    'match_data': None}
        return out_dict
        
    def matches_from_allele_mask(self, match_allele_mask) :
        series_match_masks = np.logical_and(self.allele_masks, match_allele_mask)
        return series_match_masks
    
    def series_mask_for_superset_matches(self, match_alleles_per_series, super_allele_count, min_match=0.9) :
        min_super_match_count = min_match*float(super_allele_count)
        return match_alleles_per_series >= min_super_match_count
    
    def series_mask_for_subset_matches(self, match_alleles_per_series, min_match=0.9) :
        min_sub_mask_counts = min_match*self.alleles_per_series.astype('f4')
        return match_alleles_per_series >= min_sub_mask_counts
        
    def to_out_match_data(self, series_match_masks, match_series_mask) :
        matches_size = match_series_mask.sum()
        match_data = np.zeros(matches_size, self.match_data_dtype)
        out_data = {'series_data': self.series_data[match_series_mask],
                    'sample_data': self.sample_data[match_series_mask]}
        match_data['data_index'] = out_data['series_data']['data_index']
        match_data['match_allele_mask'] = series_match_masks[match_series_mask]
        out_data['match_data'] = match_data
        return out_data
        

    def sub_super_data_from_allele_mask(self, match_allele_mask, min_match=0.9) :
        series_match_masks = self.matches_from_allele_mask(match_allele_mask)
        super_allele_count = match_allele_mask.sum()
        match_alleles_per_series = series_match_masks.sum(axis=1)
        super_match_mask = self.series_mask_for_superset_matches(match_alleles_per_series, 
                                                                 super_allele_count, min_match)
        sub_match_mask = self.series_mask_for_subset_matches(match_alleles_per_series, min_match)
        match_allele_count = match_allele_mask.sum()
        m_super = self.alleles_per_series > match_allele_count
        m_sub = np.logical_not(m_super)
        match_super = np.logical_and(m_super, super_match_mask)
        match_sub = np.logical_and(m_sub, sub_match_mask)
        match_mask = np.logical_or(match_super, match_sub)
        return self.to_out_match_data(series_match_masks, match_mask)
        
    def sub_super_data_from_data_index(self, data_index, min_match=0.9) :
        array_index = self.data_indexes.searchsorted(data_index)
        data_index_allele_mask = self.allele_masks[array_index]
        return self.sub_super_data_from_allele_mask(data_index_allele_mask, min_match)

    def hierarchy_data_from_plt_obj(self, po) :
        array_index = self.data_indexes.searchsorted(po.hierarchy_data_index)
        po.match_test_allele_mask = self.allele_masks[array_index]
        po.series_sample_match_data = self.sub_super_data_from_allele_mask(po.match_test_allele_mask, po.min_match)

    def sub_super_data_from_yes_no_indexes(self, yes_indexes=None, no_indexes=None, min_match=0.9) :
        match_allele_mask = self.yes_no_data_indexes(yes_indexes, no_indexes)
        return self.sub_super_data_from_allele_mask(match_allele_mask, min_match)

    def sub_super_yes_no_from_plt_obj(self, po) :
        po.match_test_allele_mask = self.yes_no_data_indexes(po.yes_series_data_indexes, po.no_series_data_indexes)
        po.series_sample_match_data = self.sub_super_data_from_allele_mask(po.match_test_allele_mask, po.min_match)

    def sub_super_allele_mask_yes_no_from_plt_obj(self, po) :
        if (po.yes_series_data_indexes is None) and (po.no_series_data_indexes is None) :
            po.match_test_allele_mask = po.yes_allele_mask
        else :
            po.match_test_allele_mask = self.yes_no_data_indexes(po.yes_series_data_indexes, po.no_series_data_indexes)
            po.match_test_allele_mask = np.logical_and(po.match_test_allele_mask, po.yes_allele_mask)
        po.series_sample_match_data = self.sub_super_data_from_allele_mask(po.match_test_allele_mask, po.min_match)        

    def series_start_aligned_mask(self, first_start_pos, last_start_pos) :
        m = self.series_first_pos >= first_start_pos
        m = np.logical_and(m, self.series_first_pos <= last_start_pos)
        return m

    def series_end_aligned_mask(self, first_end_pos, last_end_pos) :
        m = self.series_last_pos >= first_end_pos
        m = np.logical_and(m, self.series_last_pos <= last_end_pos)
        return m

    def series_sample_data_for_aligned_starts(self, start_pos, max_delta=5000) :
        first_start_pos = start_pos - max_delta
        last_start_pos = start_pos + max_delta
        align_mask = self.series_start_aligned_mask(first_start_pos, last_start_pos)
        return self.out_data_from_data_mask(align_mask)

    def series_sample_data_for_aligned_ends(self, end_pos, max_delta=5000) :
        first_end_pos = end_pos - max_delta
        last_end_pos = end_pos + max_delta
        align_mask = self.series_end_aligned_mask(first_end_pos, last_end_pos)
        return self.out_data_from_data_mask(align_mask)

    def or_allele_mask(self, data_indexes) :
        interval_indexes = self.data_indexes.searchsorted(data_indexes)
        out_mask = np.zeros(self.total_allele_count, '?')
        for index in interval_indexes :
            index_mask = self.allele_masks[index]
            out_mask = np.logical_or(out_mask, index_mask)
        return out_mask
        
    def and_allele_mask(self, data_indexes) :
        interval_indexes = self.data_indexes.searchsorted(data_indexes)
        out_mask = self.allele_masks[interval_indexes[0]].copy()
        for index in interval_indexes[1:] :                    
            index_mask = self.allele_masks[index]
            out_mask = np.logical_and(out_mask, index_mask)
        return out_mask

    def or_and_not_allele_mask(self, or_data_indexes, not_or_data_indexes) :
        or_allele_mask = self.or_allele_mask(or_data_indexes)
        not_or_allele_mask = self.or_allele_mask(not_or_data_indexes)
        not_or_allele_mask = np.logical_not(not_or_allele_mask)
        return np.logical_and(or_allele_mask, not_or_allele_mask)

    def and_and_not_or_allele_mask(self, and_data_indexes, not_or_data_indexes) :
        and_allele_mask = self.and_allele_mask(and_data_indexes)
        not_or_allele_mask = self.or_allele_mask(not_or_data_indexes)
        not_or_allele_mask = np.logical_not(not_or_allele_mask)
        return np.logical_and(and_allele_mask, not_or_allele_mask)

    def not_or_allele_mask(self, not_data_indexes) :
        allele_mask = self.or_allele_mask(not_data_indexes)
        return np.logical_not(allele_mask)
        
    def not_and_allele_mask(self, not_data_indexes) :
        allele_mask = self.and_allele_mask(not_data_indexes)
        return np.logical_not(allele_mask)
        
    def yes_no_data_indexes(self, yes_data_indexes=None, no_data_indexes=None) :
        if yes_data_indexes is None :
            result = self.not_or_allele_mask(no_data_indexes)
        elif no_data_indexes is None :
            result = self.and_allele_mask(yes_data_indexes)
        else :
            result = self.and_and_not_or_allele_mask(yes_data_indexes, no_data_indexes)
        result = result.astype('?')
        return result

    def yes_allele_mask_or_indexes(self, yes_allele_mask, data_indexes) :
        or_allele_mask = self.or_allele_mask(data_indexes)
        return np.logical_and(yes_allele_mask, or_allele_mask)
        
    def yes_allele_mask_and_indexes(self, yes_allele_mask, data_indexes) :
        and_allele_mask = self.and_allele_mask(data_indexes)
        return np.logical_and(yes_allele_mask, and_allele_mask)
 
    def yes_allele_mask_not_or_indexes(self, yes_allele_mask, no_data_indexes) :
        no_allele_mask = self.not_or_allele_mask(no_data_indexes)
        return np.logical_and(yes_allele_mask, no_allele_mask)
        
    def yes_allele_mask_not_and_indexes(self, yes_allele_mask, no_data_indexes) :
        no_allele_mask = self.not_and_allele_mask(no_data_indexes)
        return np.logical_and(yes_allele_mask, no_allele_mask)

    def allele_mask_from_data_index(self, data_index) :
        interval_index = self.data_indexes.searchsorted(data_index)
        out_mask = self.allele_masks[interval_index].copy()
        out_mask = out_mask.astype('?')
        return out_mask
        
    def superset_series_mask_from_allele_indexes(self, allele_indexes, min_match=0.9) :
        test_allele_masks = self.allele_masks[:,allele_indexes]
        match_alleles_per_series = test_allele_masks.sum(axis=1)
        match_thresh = min_match*float(allele_indexes.size)
        mask = match_alleles_per_series >= match_thresh
        return mask, match_alleles_per_series
        
    def subset_series_mask_from_allele_indexes(self, allele_indexes, min_match=0.9) :
        test_allele_masks = self.allele_masks[:,allele_indexes]
        match_alleles_per_series = test_allele_masks.sum(axis=1)
        mask = match_alleles_per_series >= min_match*self.alleles_per_series.astype('f4')
        return mask, match_alleles_per_series

    def superset_data_from_data_index(self, data_index, min_match=0.9) :
        index = self.data_indexes.searchsorted(data_index)
        index_allele_mask = self.allele_masks[index]
        index_allele_indexes = np.where(index_allele_mask)[0]
        data_mask, match_alleles_per_series = self.superset_series_mask_from_allele_indexes(index_allele_indexes, min_match)
        #data_mask[index] = False
        series_data = self.series_data[data_mask]
        allele_data = self.allele_masks[data_mask]
        match_alleles_per_series = match_alleles_per_series[data_mask]
        return series_data, allele_data, match_alleles_per_series
        
    def subset_data_from_allele_mask(self, allele_mask) :        
        allele_indexes = np.where(allele_mask)[0]
        data_mask, match_alleles_per_series = self.subset_series_mask_from_allele_indexes(allele_indexes)
        series_data = self.series_data[data_mask]
        allele_data = self.allele_masks[data_mask]
        match_alleles_per_series = match_alleles_per_series[data_mask]
        return series_data, allele_data, match_alleles_per_series

    def subset_data_from_data_index(self, data_index) :
        index = self.data_indexes.searchsorted(data_index)
        index_allele_mask = self.allele_masks[index]
        return self.subset_data_from_allele_mask(index_allele_mask)

                
    def superset_data_from_match_allele_mask(self, match_allele_mask, min_match=0.9) :
        match_allele_count = match_allele_mask.sum()
        if match_allele_count == 0 :
            matched_series = np.empty(0, self.series_data.dtype)
            matched_series_allele_masks = np.empty((0, self.total_allele_count), '?')
            series_match_masks = np.empty((0, self.total_allele_count), '?')
        else :
            series_match_masks  = np.logical_and(self.allele_masks, match_allele_mask)
            match_alleles_per_series = series_match_masks.sum(axis=1)
            min_match_count = min_match*float(match_allele_count)
            match_mask = match_alleles_per_series >= min_match_count
            matched_series = self.series_data[match_mask]
            matched_series_allele_masks = self.allele_masks[match_mask]
            series_match_masks = series_match_masks[match_mask]
        return matched_series, matched_series_allele_masks, series_match_masks

    def superset_data_from_yes_no_indexes(self, yes_indexes, no_indexes=None, min_match=0.9) :
        match_allele_mask = self.yes_no_data_indexes(yes_indexes, no_indexes)
        matched_series, matched_series_allele_masks, series_match_masks = (
                                      self.superset_data_from_match_allele_mask(match_allele_mask, min_match))
        return match_allele_mask, matched_series, matched_series_allele_masks, series_match_masks
        

    def subset_data_from_match_allele_mask(self, match_allele_mask, min_match=0.9) :
        series_match_masks  = np.logical_and(self.allele_masks, match_allele_mask)
        match_alleles_per_series = series_match_masks.sum(axis=1)
        min_match_counts = min_match*self.alleles_per_series.astype('f4')
        match_mask = match_alleles_per_series >= min_match_counts
        matched_series = self.series_data[match_mask]
        matched_series_allele_masks = self.allele_masks[match_mask]
        series_match_masks = series_match_masks[match_mask]
        return matched_series, matched_series_allele_masks, series_match_masks


    def subset_data_from_yes_no_indexes(self, yes_indexes, no_indexes=None, min_match=0.9) :
        match_allele_mask = self.yes_no_data_indexes(yes_indexes, no_indexes)
        return self.subset_data_from_match_allele_mask(match_allele_mask, min_match)
        




        