# -*- coding: utf-8 -*-

import tables as tb
import numpy as np
import os
from .chrom_snp_data_rdr import chrom_snp_item_factory_cls

s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)


class chrom_snp_series_tables_cls(object) :
    series_table_name_end = '_snp_series'
    series_items_table_name_end = '_snp_series_items'
    file_name_end = '_snp_series.h5'
    def __init__(self, chrom) :
        self.chrom = chrom
        table_name_start = 'chrom_' + str(chrom)
        series_table_name = table_name_start + self.series_table_name_end
        series_items_table_name = table_name_start + self.series_items_table_name_end
        file_name = 'chrom_' + str(chrom) + self.file_name_end
        file_path = os.path.join(mod_dir, file_name)
        self.h5 = tb.open_file(file_path, 'r')
        self.series_table = getattr(self.h5.root, series_table_name)
        self.series_items_table = getattr(self.h5.root, series_items_table_name)        
                
    def close(self) :
        if self.h5 is not None :
            self.h5.close()
            self.h5 = None
            self.series_table = None
            self.series_item_table = None
            
    def __del__(self) :
        self.close()

'''
class snp_item_cls(object) :
    def __init__(self, item_data, item_allele_mask) :
        self.data = item_data
        self.allele_mask = item_allele_mask
'''        
class snp_series_data_cls(object) :
    def __init__(self, series_snp_objs) :
        item_data = []
        allele_masks = []
        for obj in series_snp_objs :
            item_data.append(obj.data)
            allele_masks.append(obj.allele_mask)
        self.snp_data = np.array(item_data, item_data[0].dtype)            
        self.allele_masks = np.array(allele_masks, allele_masks[0].dtype)
        #self.allele_masks.reshape(self.snp_data.size, -1)
        
        
class chrom_snp_series_factory_cls(object) :
    def __init__(self, chrom) :
        self.chrom = chrom
        self.data_tables = chrom_snp_series_tables_cls(chrom)
        self.item_objs = chrom_snp_item_factory_cls(chrom)

    def snp_series_data(self, item_data_start, item_count) :
        item_data_bound = item_data_start + item_count
        series_items = self.data_tables.series_items_table[item_data_start:item_data_bound]
        item_snp_indexes = series_items['snp_index']
        item_objs = self.item_objs.snp_item_objs(item_snp_indexes)
        return snp_series_data_cls(item_objs)

    def item_objs_from_series_data(self, series_data) :
        item_data_start = series_data['item_data_start']
        item_count = series_data['item_count']
        return self.snp_series_data(item_data_start, item_count)
        
    def close(self) :
        if self.data_tables is not None :
            self.data_tables.close()
            self.data_tables = None
            self.item_objs.close()
            self.item_objs = None
            
    def __del__(self) :
        self.close()

