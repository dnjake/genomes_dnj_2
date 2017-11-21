# -*- coding: utf-8 -*-
from sorted_autosome_data.chrom_snp_series_data_rdr import chrom_snp_series_data_tables_cls

chrom = 2
data_tables = chrom_snp_series_data_tables_cls(chrom)

#series_data = data_tables.series_data_table[:]
#data_tables.close()



'''
>>> series_data.size
612399
>>> series_data.dtype
dtype([('data_index', '<u4'), ('first_snp_index', '<u4'), ('chrom', '<u2'), ('first_pos', '<u4'), ('last_pos', '<u4'), ('item_data_start', '<u4'), ('item_count', '<u2'), ('p90_allele_count', '<u2'), ('afr_obs_to_pred', '<f4'), ('afx_obs_to_pred', '<f4'), ('amr_obs_to_pred', '<f4'), ('eas_obs_to_pred', '<f4'), ('eur_obs_to_pred', '<f4'), ('sas_obs_to_pred', '<f4'), ('sax_obs_to_pred', '<f4')])
>>> series_first_pos = series_data['first_pos']
>>> low_pos = series_first_pos[:-1]
>>> high_pos = series_first_pos[1:]
>>> m = low_pos < high_pos
>>> m.all()
True
>>> m_single = series_data['item_count'] == 1
>>> m_single.sum()
415644
>>> import numpy as np
>>> m_multiple = np.logical_not(m_single)
>>> m_multiple.sum()
196755

'''
