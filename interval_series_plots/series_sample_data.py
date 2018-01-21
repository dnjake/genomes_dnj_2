# -*- coding: utf-8 -*-

import numpy as np
import PyQt5.QtGui as gui
from ..autosome_snp_data.allele_country_region_rdr import country_region_alleles_cls
cra = country_region_alleles_cls()
region_hues = np.array([0, 300, 60, 120, 180 ], dtype='i4')

def get_color_for_allele_mask(allele_mask) :
    allele_indexes = np.where(allele_mask)[0]
    simple_region_stats = cra.simple_region_stats(allele_indexes)
    obs_to_pred = simple_region_stats['obs_to_pred']
    ind_max_obs_to_pred = obs_to_pred.argmax()
    max_obs_to_pred = obs_to_pred[ind_max_obs_to_pred]
    if max_obs_to_pred > 2.5 :
        sat = 255
        val = 255
    elif max_obs_to_pred > 1.5 :
        sat = 128
        val = 196
    else :
        sat = 0
        val = 128
    qcv = gui.QColor.fromHsv(region_hues[ind_max_obs_to_pred], sat, val)
    hrgb = hex(qcv.rgb())
    bokeh_color = '#' + hrgb[4:10]
    return bokeh_color


def get_region_stats_from_allele_masks(allele_masks) :
    series_region_counts = cra.region_counts_from_allele_masks(allele_masks)
    series_obs_to_pred = cra.region_obs_to_pred_from_region_counts(series_region_counts)
    return series_region_counts, series_obs_to_pred

