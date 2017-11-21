# -*- coding: utf-8 -*-

'''
Two changes.
First I want all of the tables to be sorted by increasing first
snp index.  Clearly the snp indexes in the item table will not
be sorted
Second
I want all snps in the series tables because any kind of screen has
to consider single snps as a series
'''


import numpy as np
import tables as tb
#from genome_data.sample_data_plus import sample_count_cls
#from collections import namedtuple
from genomes_dnj.autosome_snp_data.chrom_snp_data_rdr import chrom_snp_data_tables_cls


'''
I think I will do only one round of matching in the first pass
and then merge series
'''

class work_snp_series_cls(object) :
    match_alleles_indexes = None
    snp_indexes = None
    data_indexes = None
    #def __init__(self, data_index, snp_index, snp_alleles_mask):

class chrom_snp_series_finder_cls(object) :
    segment_base_count = 1000000
    expr_threshold=0.9 
    specific_threshold=0.9    
    chrom_snp_pos_dtype = np.dtype([('snp_index', 'u4'), ('pos', 'u4')])
    chrom_segment_index_dtype = np.dtype([('start_snp_index', 'u4'), ('segment_snp_count', 'u4')])
    segment_snp_allele_dtype = np.dtype([('snp_index', 'u4'), ('allele_mask', 'u1', (5008,))])

    def __init__(self, chrom, data_table_writer) :
        self.chrom = chrom
        self.snp_data_rdr = chrom_snp_data_tables_cls(chrom)
        self.data_table_writer = data_table_writer
        self.segment_snp_allele_data = None
        self.segment_snp_indexes = None
        self.segment_snp_allele_masks = None
        self.last_segment_snp_series_objs = None
        self.this_segment_snp_series_objs = None
        self.last_first_snp_index = -1
        self.last_written_first_snp_index = -1
        
    # will keep track of data indexes in a segment that are still not part of a series
    # As series are extended data indexes will be removed from the set

    def find_series_members(self, snp_series_obj, match_candidate_data_indexes) :
        '''
        Why can't I do this with a vector
        There are two comparisons
        One is the per cent of the matching snp alleles that have to match
        The other is the precent of the selectors alleles that have to match
        '''
        select_allele_indexes = snp_series_obj.match_allele_indexes
        expr_threshold = self.expr_threshold*float(select_allele_indexes.size)
        specific_threshold = self.segment_alleles_specific_threshold[match_candidate_data_indexes]
        candidate_allele_masks = self.segment_snp_allele_masks[match_candidate_data_indexes]
        allele_match_masks = candidate_allele_masks[:, select_allele_indexes]
        allele_match_counts = np.sum(allele_match_masks, axis=1)
        match_mask = np.logical_and(allele_match_counts >= expr_threshold, allele_match_counts >= specific_threshold)
        match_data_indexes = match_candidate_data_indexes[match_mask]
        if match_data_indexes.size > 0 :
            snp_series_obj.data_indexes.extend(match_data_indexes.tolist())
            snp_series_obj.snp_indexes.extend(self.segment_snp_indexes[match_data_indexes].tolist())
        
    
    def find_series_candidates(self, match_allele_count, search_data_indexes) :
        '''
        This method will find the data indexes that have an allele count that could
        allow a match
        0.9*allele_count must be less then the match count
        0.9*match_count must be less then the allele_count
        '''
        match_allele_threshold = self.expr_threshold*float(match_allele_count)        
        candidate_alleles_per_snp = self.segment_alleles_per_snp[search_data_indexes]
        specific_thresholds_per_snp = self.segment_alleles_specific_threshold[search_data_indexes]
        candidate_mask = np.logical_and(candidate_alleles_per_snp >= match_allele_threshold, 
                                        specific_thresholds_per_snp <= match_allele_count)
        return search_data_indexes[candidate_mask]                                        
    
    def try_extend_last_match(self, snp_series_obj) :
        last_data_index = snp_series_obj.data_indexes[-1]
        self.search_data_indexes = np.setdiff1d(self.search_data_indexes, snp_series_obj.data_indexes)
        match_allele_mask = self.segment_snp_allele_masks[last_data_index]
        snp_series_obj.match_allele_indexes = np.where(match_allele_mask)[0]
        snp_series_obj.data_indexes = []
        match_allele_count = snp_series_obj.match_allele_indexes.size
        match_candidate_data_indexes = self.find_series_candidates(match_allele_count, self.search_data_indexes)
        self.find_series_members(snp_series_obj, match_candidate_data_indexes)
        if len(snp_series_obj.data_indexes) > 0 :
            return True
        else :
            return False

        
    def extend_old_series(self) :
        '''
        Try to extend the old snp series and write out unique snps or
        snp series for those that find no new match in the new segment
        '''
        for snp_series_obj in self.last_segment_snp_series_objs :
            #self.extend_series(obj)
            snp_series_obj.data_indexes = []
            match_allele_count = snp_series_obj.match_allele_indexes.size
            match_candidate_data_indexes = self.find_series_candidates(match_allele_count, self.search_data_indexes)
            self.find_series_members(snp_series_obj, match_candidate_data_indexes)
            if len(snp_series_obj.data_indexes) > 0 :
                while self.try_extend_last_match(snp_series_obj) :
                    continue
            out = self.data_table_writer
            series_snp_indexes = snp_series_obj.snp_indexes
            last_written_first_snp_index = snp_series_obj.snp_indexes[0]
            if last_written_first_snp_index <= self.last_written_first_snp_index :
                print 'last written', last_written_first_snp_index, self.last_written_first_snp_index
            self.last_written_first_snp_index = last_written_first_snp_index
            series_snp_indexes.sort()
            out.add_snp_series(series_snp_indexes)


    def find_all_new_series(self) :
        '''
        finds the new series for snps that were not matched by extensions to the
        old series
        '''
        while self.search_data_indexes.size > 0 :
            snp_series_obj = work_snp_series_cls()
            self.search_data_indexes.sort()
            match_data_index = self.search_data_indexes[0]
            #self.search_data_indexes = self.search_data_indexes[1:]
            match_snp_index = self.segment_snp_indexes[match_data_index]
            match_allele_mask = self.segment_snp_allele_masks[match_data_index]
            snp_series_obj.match_alleles_indexes = np.where(match_allele_mask)[0]
            snp_series_obj.data_indexes = [match_data_index]
            snp_series_obj.snp_indexes = [match_snp_index]
            if match_snp_index <= self.last_first_snp_index :
                print 'last first', match_snp_index, self.last_first_snp_index
            self.last_first_snp_index = match_snp_index
            self.this_segment_snp_series_objs.append(snp_series_obj)
            while self.try_extend_last_match(snp_series_obj) :
                continue
                         
        
    def find_segment_snp_series(self) :
        '''
        First try to extend old series.  Then look for new series
        '''
        if (self.last_segment_snp_series_objs is not None and
            (len(self.last_segment_snp_series_objs) > 0) ) :
            self.extend_old_series()        
        self.find_all_new_series()            
        
    def build_segment_data(self, start_snp_index, count) :
        '''
        Reads in the sample data for the segment and creates the search indexes
        for this segment
        '''
        bound_snp_index = start_snp_index + count
        bitpacked_allele_table = self.snp_data_rdr.snp_bitpacked_allele_values_table
        bitpacked_allele_data = bitpacked_allele_table[start_snp_index:bound_snp_index]
        segment_snp_allele_data = np.zeros(bitpacked_allele_data.size, self.segment_snp_allele_dtype)
        segment_snp_allele_masks = segment_snp_allele_data['allele_mask']
        segment_snp_allele_masks[:] = np.unpackbits(bitpacked_allele_data['bitpacked_values'], axis=1)
        segment_snp_allele_data['snp_index'] = bitpacked_allele_data['snp_index']
        self.segment_snp_allele_data = segment_snp_allele_data
        self.segment_snp_indexes = self.segment_snp_allele_data['snp_index']
        self.segment_snp_allele_masks = self.segment_snp_allele_data['allele_mask']
        self.segment_alleles_per_snp = np.sum(self.segment_snp_allele_masks, axis=1)
        self.segment_alleles_specific_threshold = self.specific_threshold*self.segment_alleles_per_snp.astype('f4')
        self.search_data_indexes = np.arange(self.segment_snp_indexes.size)

    def build_by_pos_segment_index(self) :
        '''
        Builds an index of the snps in each segment from the segment_size
        '''        
        snp_data_table = self.snp_data_rdr.snp_data_table
        #chrom_snp_data = np.zeros(snp_data_table.nrows, self.chrom_snp_pos_dtype)
        chrom_snp_data = np.zeros(snp_data_table.nrows, self.chrom_snp_pos_dtype)
        chrom_snp_data['snp_index'] = snp_data_table.col('snp_index')
        chrom_snp_data['pos'] = snp_data_table.col('pos')
        chrom_snp_pos = chrom_snp_data['pos']
        max_pos = chrom_snp_pos[-1]
        segment_starts = np.arange(0, max_pos, self.segment_base_count, dtype='u4')
        segment_start_indexes = chrom_snp_pos.searchsorted(segment_starts)
        chrom_segment_index = np.zeros(segment_starts.size, self.chrom_segment_index_dtype)
        chrom_segment_index['start_snp_index'] = segment_start_indexes
        chrom_segment_snp_counts = chrom_segment_index['segment_snp_count']
        chrom_segment_snp_counts[:-1] = segment_start_indexes[1:] - segment_start_indexes[:-1]
        chrom_segment_snp_counts[-1] = chrom_snp_data.size - segment_start_indexes[-1]
        self.chrom_segment_index = chrom_segment_index

    def initialize_segment(self, segment_start_snp_index, segment_snp_count) :
        self.build_segment_data(segment_start_snp_index, segment_snp_count)
        self.last_segment_snp_series_objs = self.this_segment_snp_series_objs
        self.this_segment_snp_series_objs = []

    def write_this_segment_objs(self) :
        '''
        Handle the objs from the last segment
        '''
        for snp_series_obj in self.this_segment_snp_series_objs :
            out = self.data_table_writer
            series_snp_indexes = snp_series_obj.snp_indexes
            series_snp_indexes.sort()
            out.add_snp_series(series_snp_indexes)

    def organize_chrom_processing(self) :
        '''
        handles first, last, and segment iteration
        '''
        self.build_by_pos_segment_index()
        self.this_segment_snp_series_objs = []
        self.last_segment_snp_series_objs = []
        for segment_index_item in self.chrom_segment_index :
            start_snp_index, count = segment_index_item
            print start_snp_index, count
            if count == 0 :
                self.write_this_segment_objs()
                self.this_segment_snp_series_objs = []
                continue
            self.initialize_segment(start_snp_index, count)
            self.find_segment_snp_series()
        self.write_this_segment_objs()
        self.data_table_writer.close()
        self.snp_data_rdr.close()                

class snp_series_table_writer_cls(object) :
    filters = tb.Filters(complevel=5, complib='zlib')
    file_name_end = '_snp_series.h5'
    series_table_name_end = '_snp_series'
    series_item_table_name_end = '_snp_series_items'
    series_table_dtype = np.dtype([('first_snp_index', 'u4'), ('item_data_start', 'u4'), ('item_count', 'u2') ])
    series_item_table_dtype = np.dtype([('snp_index', 'u4'), ('first_snp_index', 'u4')])
    series_table_buffer_size = 50000
    series_item_table_buffer_size = 200000
    
    def __init__(self, chrom) :
        table_name_start = 'chrom_' + str(chrom)
        file_path = table_name_start + self.file_name_end
        self.series_table_name = table_name_start + self.series_table_name_end
        self.series_item_table_name = table_name_start + self.series_item_table_name_end
        self.h5 = tb.open_file(file_path, 'w', filters=self.filters)
        self.series_table = self.h5.create_table('/', self.series_table_name, description=self.series_table_dtype)
        self.series_item_table = self.h5.create_table('/', self.series_item_table_name, description=self.series_item_table_dtype)        
        self.series_table_buffer = np.zeros(self.series_table_buffer_size, self.series_table_dtype)
        self.series_item_table_buffer = np.zeros(self.series_item_table_buffer_size, self.series_item_table_dtype)
        self.series_item_table_nrows = 0
        self.series_table_top = 0
        self.series_item_table_top = 0
        
    def write_series_table(self, replace=True) :
        if self.series_table_top > 0 :
            self.series_table.append(self.series_table_buffer[:self.series_table_top])
        if replace :
            self.series_table_buffer = np.zeros(self.series_table_buffer_size, self.series_table_dtype)
            self.series_table_top = 0
        else :
            self.series_table.close()
            self.series_table_buffer = None
            self.series_table = None
        
    def write_series_item_table(self, replace=True) :
        if self.series_item_table_top > 0 :
            self.series_item_table.append(self.series_item_table_buffer[:self.series_item_table_top])
        if replace :
            self.series_item_table_buffer = np.zeros(self.series_item_table_buffer_size, self.series_item_table_dtype)
            self.series_item_table_top = 0
        else :
            self.series_item_table.close()
            self.series_item_table_buffer = None
            self.series_item_table = None

            
    def add_snp_series(self, series_snp_indexes) :
        first_snp_index = series_snp_indexes[0]
        series_length = len(series_snp_indexes)
        snp_series = (first_snp_index, self.series_item_table_nrows, series_length)
        self.series_table_buffer[self.series_table_top] = snp_series
        self.series_table_top += 1
        if self.series_table_top == self.series_table_buffer_size :
            self.write_series_table()
        after_items_top = self.series_item_table_top + series_length
        if after_items_top > self.series_item_table_buffer_size :
            self.write_series_item_table()        
        for snp_index in series_snp_indexes :
            self.series_item_table_buffer[self.series_item_table_top] = (snp_index, first_snp_index)
            self.series_item_table_top += 1
        self.series_item_table_nrows += series_length

    def close(self) :
        if self.h5 is not None :
            self.write_series_table(False)
            self.write_series_item_table(False)
            self.h5.close()
            self.h5 = None

    def __del__(self) :
        if self.h5 is not None :
            self.close()













































