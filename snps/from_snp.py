# -*- coding: utf-8 -*-
import numpy as np
from ..autosome_snp_data import chrom_snp_data_rdr as sdr
from ..autosome_snp_data import chrom_snp_series_rdr as ssr
from ..autosome_snp_data import chrom_snp_series_data_rdr as ssdr

class from_snp_data_cls(object) :
    def __init__(self, chrom) :
        self.chrom = chrom
        self.snp_rdr_obj = sdr.chrom_snp_data_tables_cls(self.chrom)
        self.snp_table = self.snp_rdr_obj.snp_data_table
        
    def snp_data_from_id(self, snp_id) :
        rows = self.snp_table.where('id == snp_id')
        try :
            row = rows.next()
        except StopIteration :
            return None
        return row['snp_index'], row['pos'], row['id']
        
        
    def close(self) :
        if self.snp_rdr_obj is not None :
            self.snp_rdr_obj.close()
            self.snp_rdr_obj = None
            self.snp_table = None

    def __del__(self) :
        self.close()
        
        
class from_snp_series_cls(object) :
    def __init__(self, chrom) :
        self.chrom = chrom
        self.snp_series_rdr = ssr.chrom_snp_series_tables_cls(chrom)
        self.series_items_table = self.snp_series_rdr.series_items_table
        
    def series_first_index_from_snp_index(self, index) :
        rows = self.series_items_table.where('snp_index == index')
        try :
            row = rows.next()
        except StopIteration :
            return None
        return row['first_snp_index']
        
    def close(self) :
        if self.snp_series_rdr is not None :
            self.snp_series_rdr.close()
            self.snp_series_rdr = None
            self.series_items_table = None

    def __del__(self) :
        self.close()
        
class from_snp_series_data_cls(object) :
    def __init__(self, chrom) :
        self.chrom = chrom
        self.series_data_rdr = ssdr.chrom_snp_series_data_tables_cls(chrom)
        self.data_table = self.series_data_rdr.series_data_table
        
    def series_data_from_first_snp_index(self, index) :
        rows = self.data_table.where('first_snp_index == index')
        try :
            row = rows.next()
        except StopIteration :
            return None
        return row[:]

    def close(self) :
        if self.series_data_rdr is not None :
            self.series_data_rdr.close()
            self.series_data_rdr = None
            self.data_table = None

    def __del__(self) :
        self.close()

def series_data_from_ids(chrom, snp_ids) :
    fsdo = from_snp_data_cls(chrom)
    fsso = from_snp_series_cls(chrom)
    fssdo = from_snp_series_data_cls(chrom)
    out_data = []
    for snp_id in snp_ids :
        snp_data = fsdo.snp_data_from_id(snp_id)
        if snp_data is None :
            out_data.append((snp_id, None, None))
        else :
            snp_index = snp_data[0]
            first_snp_index = fsso.series_first_index_from_snp_index(snp_index)
            series_data = fssdo.series_data_from_first_snp_index(first_snp_index)
            series_index = series_data[0]
            out_data.append((snp_id, snp_index, series_index))
    fsdo.close()
    fsso.close()
    fssdo.close()
    return out_data

def print_series_from_ids(chrom, snp_ids) :
    snp_series_data = series_data_from_ids(chrom, snp_ids)
    print snp_series_data

snp_and_series_data_dtype = np.dtype([('snp_index', 'u4'), ('id', 'S11'), ('pos', 'u4'), ('series_index', 'u4'),
                                      ('first_pos', 'u4'), ('length', 'u4'), ('snp_count', 'u2'), ('allele_count', 'u2')])
    
def snp_and_series_data_from_indexes(chrom, snp_items) :
    snp_rdr_obj = sdr.chrom_snp_data_tables_cls(chrom)
    snp_table = snp_rdr_obj.snp_data_table
    series_data_rdr = ssdr.chrom_snp_series_data_tables_cls(chrom)
    series_table = series_data_rdr.series_data_table
    out_data = []
    for snp_id, snp_index, series_index in snp_items :
        snp_data = snp_table[snp_index]
        id = snp_data['id']
        pos = snp_data['pos']
        series_data = series_table[series_index]
        first_pos = series_data['first_pos']
        last_pos = series_data['last_pos']
        length = last_pos - first_pos
        snp_count = series_data['item_count']
        allele_count = series_data['p90_allele_count']
        out_data.append((snp_index, id, pos, series_index, first_pos, length, snp_count, allele_count))
    snp_rdr_obj.close()
    series_data_rdr.close()
    out_data = np.array(out_data, dtype=snp_and_series_data_dtype)
    return out_data        
    

        