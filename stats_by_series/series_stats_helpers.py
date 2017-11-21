# -*- coding: utf-8 -*-

import numpy as np




def lengths_mean_median_max(series_data) :
    series_lengths = series_data['last_pos'] - series_data['first_pos']
    mean_series_length = np.mean(series_lengths)
    median_series_length = np.median(series_lengths)
    max_series_length = series_lengths.max()
    return mean_series_length, median_series_length, max_series_length

def samples_mean_median(series_data) :
    series_samples = series_data['p90_allele_count']
    mean_samples_per_series = np.mean(series_samples)
    median_samples_per_series = np.median(series_samples)
    return mean_samples_per_series, median_samples_per_series

def snps_count_mean_median_max(series_data) :
    series_snps = series_data['item_count']
    snp_count = series_snps.sum()
    mean_snps_per_series = np.mean(series_snps)
    median_snps_per_series = np.median(series_snps)
    max_snps_per_series = series_snps.max()
    return snp_count, mean_snps_per_series, median_snps_per_series, max_snps_per_series

def series_data_no_low_some_high_all_afr(series_data) :
    data = []
    sd = series_data
    m_work = sd['afr_obs_to_pred'] < 0.1
    data.append(sd[m_work])
    m_work = sd['afr_obs_to_pred'] >= 0.1
    m_work = np.logical_and(m_work, sd['afr_obs_to_pred'] < 0.5)
    data.append(sd[m_work])
    m_work = sd['afr_obs_to_pred'] >= 0.5
    m_work = np.logical_and(m_work, sd['afr_obs_to_pred'] < 2.0)
    data.append(sd[m_work])
    m_high = sd['afr_obs_to_pred'] >= 2.0
    m_work = sd['eas_obs_to_pred'] < 0.01
    m_work = np.logical_and(m_work, sd['eur_obs_to_pred'] < 0.01)
    m_work = np.logical_and(m_work, sd['sas_obs_to_pred'] < 0.01)
    m_work = np.logical_and(m_work, sd['sax_obs_to_pred'] < 0.01)
    m_not_work = np.logical_not(m_work)
    m_high = np.logical_and(m_high, m_not_work)
    data.append(sd[m_high])
    data.append(sd[m_work])
    return data
    
def series_counts_by_pop_type(series_data) :
    pop_type_series_data = series_data_no_low_some_high_all_afr(series_data)
    series_counts = []
    for sd in pop_type_series_data :
        series_counts.append(sd.size)
    return series_counts
    
def pos_min_max_length_mean(series_data):
    pos_min = series_data['first_pos'].min()
    pos_max = series_data['last_pos'].max()    
    series_lengths = series_data['last_pos'] - series_data['first_pos']
    length_mean = series_lengths.mean()
    return pos_min, pos_max, length_mean

def samples_in_series(series_allele_masks) :
    series_per_sample = series_allele_masks.sum(axis=0)
    m = series_per_sample > 0
    return m.sum()

def series_count_sample_weighted_length(series_data) :
    series_count = series_data.size
    series_lengths = series_data['last_pos'] - series_data['first_pos']
    sample_weighted_lengths = series_data['p90_allele_count']*series_lengths
    sample_weighted_length = sample_weighted_lengths.sum()
    return series_count, sample_weighted_length

def snps_count_mean_sample_weighted(series_data) :
    snps_count = series_data['item_count'].sum()
    snps_mean = series_data['item_count'].mean()
    snps_sample_weighted = series_data['item_count']*series_data['p90_allele_count']
    snps_sample_weighted = snps_sample_weighted.sum()
    return snps_count, snps_mean, snps_sample_weighted


def mask_no_afr(series_data) :
    sd = series_data
    m = sd['afr_obs_to_pred'] < 0.1
    return m
    
def mask_low_afr(series_data) :
    sd = series_data
    m = sd['afr_obs_to_pred'] >= 0.1
    m = np.logical_and(m, sd['afr_obs_to_pred'] < 0.5)
    return m

def mask_some_afr(series_data) :
    sd = series_data
    m = sd['afr_obs_to_pred'] >= 0.5
    m = np.logical_and(m, sd['afr_obs_to_pred'] < 2.0)
    return m

def mask_high_afr(series_data, m_all_afr) :
    sd = series_data
    m = sd['afr_obs_to_pred'] >= 2.0
    m_not_all_afr =  np.logical_not(m_all_afr)
    m = np.logical_and(m, m_not_all_afr)
    return m

def mask_all_afr(series_data) :
    sd = series_data
    m = sd['eas_obs_to_pred'] < 0.01
    m = np.logical_and(m, sd['eur_obs_to_pred'] < 0.01)
    m = np.logical_and(m, sd['sas_obs_to_pred'] < 0.01)
    m = np.logical_and(m, sd['sax_obs_to_pred'] < 0.01)
    return m

def data_bin_counts(data, bins) :
    data = data.copy()
    data.sort()
    ind_bins = data.searchsorted(bins)
    all_inds = [0]
    all_inds.extend(ind_bins)
    all_inds.append(data.size)
    all_inds = np.array(all_inds, 'i4')
    bin_counts = all_inds[1:] - all_inds[:-1]
    return bin_counts


















