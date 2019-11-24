# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import numpy as np
import tables as tb

import genomes_dnj_2.autosome_snp_data.chrom_snp_data_rdr as sdr
import genomes_dnj_2.autosome_snp_data.chrom_snp_series_rdr as ssr
import genomes_dnj_2.autosome_snp_data.chrom_snp_series_data_rdr as ssdr
import genomes_dnj_2.repeat_analysis_data.alu_snp_data as asd


import genomes_dnj_2.repeat_analysis.num_16_decode as n16d
filters = tb.Filters(complevel=5, complib='zlib')


class data_writer_cls(object) :
    alu_sequence_data_folder = 'grch37_hg19_alu_data'
    repeat_snp_dtype = np.dtype([('chrom', np.uint16), ('alu_offset', np.uint16), ('snp_index', '<u4'), ('pos', '<u4'),
                                 ('repeat_index', np.uint32), ('series_index', np.uint32), ('all_count', '<u2'),
                                 ('item_count', '<u2'), ('strand', 'S1'), ('ref', 'S1'), ('alt', 'S1'),
                                 ('not_expressed_is_variant', 'u1')])
    
    def __init__(self, repeat_snp_obj) :
        self.repeat_snp_obj = repeat_snp_obj
        self.snp_data = self.repeat_snp_obj.snp_data
        self.snp_data_field_names = self.snp_data.dtype.names
        self.repeats = self.repeat_snp_obj.repeats
        
    def chrom_repeat_snp_data(self, chrom) :
        bound = chrom + 1
        inds = self.snp_data['chrom'].searchsorted([chrom, bound])
        data = self.snp_data[inds[0]:inds[1]]
        return data

    def chrom_snp_data(self, chrom, snp_indexes) :
        rdr = sdr.chrom_snp_data_tables_cls(chrom)
        snps = rdr.snp_data_table[snp_indexes]
        rdr.close()
        return snps
    
    def chrom_series_snp_first_indexes(self, chrom, snp_indexes) :
        sst = ssr.chrom_snp_series_tables_cls(chrom)
        snp_series_items = sst.series_items_table[:]
        snp_series_items.sort(order='snp_index')
        item_indexes = snp_series_items['snp_index'].searchsorted(snp_indexes)
        snp_series_items = snp_series_items[item_indexes]
        series = sst.series_table[:]
        series.sort(order='first_snp_index')
        sst.close()
        inds_snp_items = series['first_snp_index'].searchsorted(snp_series_items['first_snp_index'])
        item_first_snps = series[inds_snp_items]
        return item_first_snps
    
    def chrom_snp_series_indexes(self, chrom, series_snp_first_indexes) :
        ssdo =ssdr.chrom_snp_series_data_tables_cls(chrom)
        series_data = ssdo.series_data_table[:]
        series_data.sort(order='first_snp_index')
        inds = series_data['first_snp_index'].searchsorted(series_snp_first_indexes)
        return series_data['data_index'][inds]
    
    def do_chrom(self, chrom) :
        crsd = self.chrom_repeat_snp_data(chrom)
        csd = self.chrom_snp_data(chrom, crsd['snp_index'])
        csfi = self.chrom_series_snp_first_indexes(chrom, crsd['snp_index'])
        cd = np.zeros(crsd.size, dtype=self.repeat_snp_dtype)
        for name in self.snp_data_field_names :
            cd[name] = crsd[name]
        cd['ref'] = csd['ref']
        cd['alt'] = csd['alt']
        cd['all_count'] = csd['all_count']
        cd['not_expressed_is_variant'] = csd['not_expressed_is_variant']
        cd['item_count'] = csfi['item_count']
        cd['series_index'] = self.chrom_snp_series_indexes(chrom, csfi['first_snp_index'])
        cd['strand'] = self.strand_polarity        
        return cd
        
    def do_chroms(self) :
        out_data = []
        for chrom in range(1, 23) :
            cd = self.do_chrom(chrom)
            out_data.append(cd)
        self.repeat_full_snp_data = np.concatenate(out_data)
        
    def write_data(self) :
        file_path = os.path.join(self.alu_sequence_data_folder, self.file_name)
        h5 = tb.open_file(file_path, 'w', filters=filters)
        snp_data_table = h5.create_table('/', self.snp_data_table_name, description=self.repeat_full_snp_data.dtype)
        snp_data_table.append(self.repeat_full_snp_data)
        repeat_table = h5.create_table('/', self.alu_repeat_table_name, description=self.repeats.dtype)
        repeat_table.append(self.repeats)
        h5.close()
        
    def do_write(self) :
        self.do_chroms()
        self.write_data()
        
class pos_strand_data_writer_cls(data_writer_cls) :
    file_name = 'pos_alu_repeat_full_snps.h5'
    snp_data_table_name = 'pos_alu_repeat_full_snps'
    alu_repeat_table_name = 'alu_plus_strand_repeats'
    strand_polarity = '+'
    
    def __init__(self) :
        pos_snp_obj = asd.pos_alu_snp_data_cls()
        data_writer_cls.__init__(self, pos_snp_obj)
        

class neg_strand_data_writer_cls(data_writer_cls) :
    file_name = 'neg_alu_repeat_full_snps.h5'
    snp_data_table_name = 'neg_alu_repeat_full_snps'
    alu_repeat_table_name = 'alu_minus_strand_repeats'
    strand_polarity = '-'
    
    def __init__(self) :
        neg_snp_obj = asd.neg_alu_snp_data_cls()
        data_writer_cls.__init__(self, neg_snp_obj)
        
    
'''
wo = neg_strand_data_writer_cls()
wo.do_write()
'''
