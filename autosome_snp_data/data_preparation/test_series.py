# -*- coding: utf-8 -*-
import numpy as np

from autosome_snp_data.chrom_snp_data_rdr import chrom_snp_data_tables_cls
chrom = 2
data_tables = chrom_snp_data_tables_cls(chrom)
dt = data_tables.snp_data_table


snp_indexes = dt.col('snp_index')
pos_values = dt.col('pos')
data_indexes = np.arange(dt.nrows, dtype='u4')

data_indexes = np.arange(dt.nrows, dtype='u4')
m = snp_indexes == data_indexes
assert m.all()

pos_low = pos_values[:-1]
pos_high = pos_values[1:]
m = pos_low < pos_high
assert m.all()

bp_snp_indexes = data_tables.snp_bitpacked_allele_values_table.col('snp_index')
m = snp_indexes == bp_snp_indexes

'''
>>> m.all()
True
'''
data_tables.close()


from sorted_autosome_data.chrom_snp_series_rdr import chrom_snp_series_tables_cls
series_tables = chrom_snp_series_tables_cls(chrom)

'''
>>> series_tables.series_items_table.nrows
1547546
>>> pos_values.size
1547546

'''

items_table = series_tables.series_items_table
series_items = items_table[:]
items_first_snp_index = series_items['first_snp_index']
first_low = items_first_snp_index[:-1]
first_high = items_first_snp_index[1:]
m = first_low <= first_high

'''
>>> m.all()
False
'''
mn = np.logical_not(m)
out_of_order = first_low[mn]

'''
>>> out_of_order.size
240
'''

series_tables.close()
 
'''
This is the new data
The change was never to extend a series beyond a second segment
Most likely the change will not lose any series members.
>>> series_tables.series_items_table.nrows
1547546
>>> series_tables.series_table.nrows
612399
>>> snp_indexes.size
1547546
>>> series_tables.series_table.dtype
dtype([('first_snp_index', '<u4'), ('item_data_start', '<u4'), ('item_count', '<u2')])
>>> series_data = series_tables.series_table[:]
>>> m = series_data['item_count'] == 1
>>> m.sum()
415644
>>> ms = series_data['item_count'] > 1
>>> ms.sum()
196755
>>> series_tables.close()
>>> mn.sum()
0

'''















   