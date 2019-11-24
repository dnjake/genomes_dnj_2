# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

import os
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')

import genomes_dnj_2.repeat_analysis_data.alu_16_base_data as abd

class center_num_16_cls(object) :
        
    def __init__(self, center_num_16, num_data, num_counts) :
        self.center_num_16 = center_num_16
        self.num_data = num_data
        self.num_counts = num_counts
        self.center_inds = np.where(self.num_data['num_16'] == self.center_num_16)[0]
        repeat_indexes = self.num_data['repeat_index'][self.center_inds]
        self.unique_repeat_indexes = np.unique(repeat_indexes)
        self.center_duplicate = False
        self.center_count = self.center_inds.size
        self.min_shared_count = self.center_count/2
        
    def validate_num_16(self, exp_num_16) :
        ind = self.num_counts['num_16'].searchsorted(exp_num_16)
        num_16_item = self.num_counts[ind]
        num, count, already_explored = num_16_item
        assert num == exp_num_16
        if num == self.center_num_16 :
            self.center_duplicate = True
            return False
        if already_explored :
            return False
        num_16_item['already_explored'] = True
        return True
        
        
    def validate_offset(self, offset_num_16, offset_num_16_count, offset) :
        if offset_num_16_count >= self.min_shared_count :
            if self.validate_num_16(offset_num_16) :
                self.links.append((self.center_num_16, offset_num_16, offset, offset_num_16_count))
                return True
        return False
        
    def exp_offset_from_center(self, offset) :
        data_size = self.num_data.size
        data_offset = self.num_data['alu_offset']
        if (offset > 0) :
            offset_inds = self.center_inds + offset
            valid_mask = offset_inds < data_size
            valid_center_inds = self.center_inds[valid_mask]
            offset_inds = offset_inds[valid_mask]
        else :
            valid_mask = self.center_inds >= -offset
            valid_center_inds = self.center_inds[valid_mask]
            offset_inds = valid_center_inds + offset
        
        offset_alu_offset = data_offset[valid_center_inds] + offset
        valid_mask = data_offset[offset_inds] == offset_alu_offset        
        offset_inds = offset_inds[valid_mask]
        offset_num_16 = self.num_data['num_16'][offset_inds]
        if offset_num_16.size == 0 :
            return False
        unique_offset_num_16, counts = np.unique(offset_num_16, return_counts=True)
        arg_max_count = np.argmax(counts)
        max_offset_num_16 = unique_offset_num_16[arg_max_count]
        max_offset_num_16_count = counts[arg_max_count]
        return self.validate_offset(max_offset_num_16, max_offset_num_16_count, offset)
    
    def find_positive_links(self) :
        offset = 0
        Found = True
        while Found :
            offset += 1
            Found = self.exp_offset_from_center(offset)
            
    def find_negative_links(self) :
        offset = 0
        Found = True
        while Found :
            offset -= 1
            Found = self.exp_offset_from_center(offset)
        
    def find_links(self) :
        ind = self.num_counts['num_16'].searchsorted(self.center_num_16)
        center_item = self.num_counts[ind]
        num_16, count, already_explored = center_item
        assert num_16 == self.center_num_16
        assert self.center_count == count
        assert already_explored == False
        center_item['already_explored'] = True
        self.links = []
        self.links.append((self.center_num_16, self.center_num_16, 0, self.center_count))
        self.find_positive_links()
        self.find_negative_links()
        
class cluster_finder_cls(object) :
    alu_sequence_data_folder = 'grch37_hg19_alu_data'
    num_count_dtype = np.dtype([('num_16', np.uint32), ('total_count', np.uint32), ('already_explored', np.bool)])
    link_dtype = np.dtype([('center', np.uint32), ('linked', np.uint32), ('offset', np.int32), ('shared_repeat_count', np.uint32)])
    poly_dtype = np.dtype([('num_16', np.uint32), ('total_count', np.uint32), ('unique_repeat_count', np.uint32)])
    min_num_16_count = 100

    
    def __init__(self, data_obj) :
        self.data_obj = data_obj
        self.num_data = data_obj.num_data
        self.data_obj.read_num_counts()
        nc = self.data_obj.num_counts
        self.num_counts = np.zeros(nc.size, dtype=self.num_count_dtype)
        self.num_counts['num_16'] = nc['num_16']
        self.num_counts['total_count'] = nc['total_count']
    

    def find_clusters(self) :
        count_args = self.num_counts['total_count'].argsort()
        count_args = count_args[::-1]
        m = self.num_counts['total_count'][count_args] >= self.min_num_16_count
        count_args = count_args[m]
        pco = 0 
        for arg in count_args :
            num_16, count, already_explored = self.num_counts[arg]
            if already_explored :
                continue
            if pco % 1000 == 0 :
                print(num_16, count)
            pco += 1
            num_obj = center_num_16_cls(num_16, self.num_data, self.num_counts)
            num_obj.find_links()
            if num_obj.center_duplicate :
                data = num_obj.center_num_16, num_obj.center_count, num_obj.unique_repeat_indexes.size
                self.poly_table.append([data])
            else :
                self.clusters_table.append(num_obj.links)
    
    def open_output(self) :
        #cluster_file_path = os.path.join(self.alu_sequence_data_folder, self.cluster_file_name)
        cluster_file_path = self.cluster_file_name
        self.h5 = tb.open_file(cluster_file_path, 'w', filters=filters)
        self.clusters_table = self.h5.create_table('/', self.cluster_table_name, description=self.link_dtype)
        self.poly_table = self.h5.create_table('/', self.poly_table_name, description=self.poly_dtype)
        
    def close(self) :
        self.h5.close()
        
    def do_write(self) :
        self.open_output()
        self.find_clusters()
        self.close()

class pos_cluster_finder_cls(cluster_finder_cls) :
    cluster_file_name = 'pos_alu_num_16_clusters.h5'
    cluster_table_name = 'pos_alu_num_16_clusters'
    poly_table_name = 'pos_alu_repeat_poly_num_16s'
    
    def __init__(self) :
        data_obj = abd.pos_alu_data_cls()
        cluster_finder_cls.__init__(self, data_obj)

class neg_cluster_finder_cls(cluster_finder_cls) :
    cluster_file_name = 'neg_alu_num_16_clusters.h5'
    cluster_table_name = 'neg_alu_num_16_clusters'
    poly_table_name = 'neg_alu_repeat_poly_num_16s'
    
    def __init__(self) :
        data_obj = abd.neg_alu_data_cls()
        cluster_finder_cls.__init__(self, data_obj)


'''
cfo = neg_cluster_finder_cls()
cfo.do_write()
'''
















