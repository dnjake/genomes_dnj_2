# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import genomes_dnj_2.repeat_analysis_data.cpg_island_data as cid
import genomes_dnj_2.repeat_analysis_data.all_repeat_data as ard


class alu_cpg_finder_cls(object) :
    data_dtype = np.dtype([('repeat_index', np.uint32), ('cpg_index', np.uint32), ('repeat_name', 'S20')])
    count_dtype = np.dtype([('repeat_name', 'S20'), ('start', np.uint32), ('count', np.uint32), 
                            ('name_count', np.uint32)])    
    
    def __init__(self) :
        self.read_data()
        self.process_chroms()
        self.count_repeat_names()
        
    def read_data(self) :
        all_cpg = cid.all_cpg_island_cls()
        self.cpg = all_cpg.all_cpg_islands[1:]
        rpd = ard.alu_repeat_cls()
        self.alu_repeat_obj = rpd
        self.alu_repeats = rpd.repeat_data
        
    def chrom_cpgs(self, chrom) :
        bound = chrom + 1
        cis = self.cpg['chrom'].searchsorted([chrom, bound])
        return self.cpg[cis[0]:cis[1]]
    
    def chrom_alus(self, chrom) :
        bound = chrom + 1
        cis = self.alu_repeats['chrom'].searchsorted([chrom, bound])
        return self.alu_repeats[cis[0]:cis[1]]
    
    def process_chroms(self) :
        alu_cpg_data = []
        for chrom in range(1, 23) :
            cpgs = self.chrom_cpgs(chrom)
            alus = self.chrom_alus(chrom)
            inds_start_start = alus['start_pos'].searchsorted(cpgs['start']-1)
            inds_end_end = alus['end_pos'].searchsorted(cpgs['end'])
            m = inds_start_start > inds_end_end
            inds = inds_end_end[m]
            ind_alus = alus[inds]
            ind_cpgs = cpgs[m]
            chrom_data = np.zeros(ind_alus.size, dtype=self.data_dtype)
            chrom_data['repeat_index'] = ind_alus['index']
            chrom_data['cpg_index'] = ind_cpgs['index']
            chrom_data['repeat_name'] = ind_alus['repeat_name']
            alu_cpg_data.append(chrom_data)
        self.alu_cpg_data = np.concatenate(alu_cpg_data)
            
    def repeats_from_strand(self, strand) :
        ris = self.alu_cpg_data['repeat_index']
        iris = self.alu_repeats['index'].searchsorted(ris)
        cpg_repeats = self.alu_repeats[iris]
        m = cpg_repeats['strand'] == strand
        return cpg_repeats[m]
        
    def count_repeat_names(self) :
        self.alu_cpg_data.sort(order=['repeat_name', 'repeat_index'])
        names, starts, counts = np.unique(self.alu_cpg_data['repeat_name'], return_index=True, return_counts=True)
        cd = np.zeros(names.size, dtype=self.count_dtype)
        cd['repeat_name'] = names
        cd['start'] = starts
        cd['count'] = counts
        name_counts = cd['name_count']
        self.repeat_name_counts = cd
        for i in range(names.size) :
            name = names[i]
            name_repeats = self.alu_repeat_obj.repeats_from_name(name)
            name_counts[i] = name_repeats.size

    def repeats_from_name(self, name, strand=None) :
        ind_name = self.repeat_name_counts['repeat_name'].searchsorted(name)
        repeat_name, start, count, name_count = self.repeat_name_counts[ind_name]
        bound = start + count
        name_repeat_indexes = self.alu_cpg_data['repeat_index'][start:bound]
        name_repeats = self.alu_repeat_obj.repeats_from_indexes(name_repeat_indexes)
        if strand is not None :
            m = name_repeats['strand'] == strand
            name_repeats = name_repeats[m]
        return name_repeats
