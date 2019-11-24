# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import numpy as np
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d

comp_num_count_dtype = np.dtype([('num_16', np.uint32), ('letters', 'S16'), ('repeat_num_count', np.uint32),
                                 ('data_num_count', np.uint32), ('rel_expr', np.float32)])

alu_low_range = 1, 121
alu_high_range = 135, 301

def low_range_num_counts(num_data_obj) :
    offset_start, offset_bound = alu_low_range
    num_counts = num_data_obj.num_counts_from_offset_range(offset_start, offset_bound)
    return num_counts

def high_range_num_counts(num_data_obj) :
    offset_start, offset_bound = alu_high_range
    num_counts = num_data_obj.num_counts_from_offset_range(offset_start, offset_bound)
    return num_counts

def low_range_repeat_indexes_from_num(num_16, num_data_obj) :
    offset_start, offset_bound = alu_low_range    
    repeat_indexes = num_data_obj.repeat_indexes_from_num_and_offset_range(num_16, offset_start, offset_bound)
    return repeat_indexes

def high_range_repeat_indexes_from_num(num_16, num_data_obj) :
    offset_start, offset_bound = alu_high_range    
    repeat_indexes = num_data_obj.repeat_indexes_from_num_and_offset_range(num_16, offset_start, offset_bound)
    return repeat_indexes

def alu_name_counts_from_repeat_indexes(repeat_indexes, data_obj) :
    inds = data_obj.repeats['index'].searchsorted(repeat_indexes)
    alu_names = data_obj.repeats['repeat_name'][inds]
    names, counts = np.unique(alu_names, return_counts=True)
    return zip(names, counts)
    
def merge_num_counts(subset_num_counts, total_num_counts, subset_repeat_count, total_repeat_count, filter_100=True) :
    dnc = total_num_counts
    rnc = subset_num_counts
    da = np.zeros(dnc.size, dtype=comp_num_count_dtype)
    da['num_16'] = dnc['num_16']
    da['data_num_count'] = dnc['count']
    inds = da['num_16'].searchsorted(rnc['num_16'])
    da['repeat_num_count'][inds] = rnc['count']
    pred_factor = float(subset_repeat_count)/float(total_repeat_count)
    pred_num_counts = pred_factor*da['data_num_count'].astype(np.float32)
    da['rel_expr'] = da['repeat_num_count'].astype(np.float32)/pred_num_counts
    if filter_100 :
        m = da['data_num_count'] >= 100
        da = da[m]
    da['letters'] = n16d.num_16_strs(da['num_16'])
    return da

def counts_for_offset(parent_ndo, subset_ndo, offset) :
    nc_total = parent_ndo.num_counts_from_offset(offset)
    total_repeat_count = nc_total['count'].sum()
    nc_subset = subset_ndo.num_counts_from_offset(offset)
    subset_repeat_count = nc_subset['count'].sum()
    da = merge_num_counts(nc_subset, nc_total, subset_repeat_count, total_repeat_count)
    da.sort(order=['data_num_count', 'repeat_num_count'])
    da = da[::-1]
    return da
    
def repeats_from_indexes(repeat_indexes, repeats) :
    inds = repeats['index'].searchsorted(repeat_indexes)
    m = repeats['index'][inds] == repeat_indexes
    assert m.all()
    return repeats[inds]
    
def repeats_from_cg_counts(repeats, min_count=None, max_count=None, equal_count=None) :
    cg_counts = repeats['cg_count']
    if min_count is not None :
        m = cg_counts >= min_count
    else :
        m = cg_counts >= 0
    if max_count is not None :
        m = np.logical_and(m, cg_counts <= max_count)
    if equal_count is not None :
        m = np.logical_and(m, cg_counts == equal_count)
    return repeats[m]

cg_count_dtype = np.dtype([('cg_count', np.uint16), ('count', np.uint32)])        
def cg_counts_from_repeat_indexes(repeat_indexes, all_repeats) :
    repeats = repeats_from_indexes(repeat_indexes, all_repeats)
    cg_counts, counts = np.unique(repeats['cg_count'], return_counts=True)
    ca = np.zeros(cg_counts.size, dtype=cg_count_dtype)
    ca['cg_count'] = cg_counts
    ca['count'] = counts
    return ca
    
    
    
    
    
    
    
    
    
    
    
    
    
    