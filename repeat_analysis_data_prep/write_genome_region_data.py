# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
#import genomes_dnj_2.repeat_analysis_data.gene_transcript_data as gtd
#import genomes_dnj_2.repeat_analysis_data.all_repeat_data as ard
#import genomes_dnj_2.repeat_analysis_data.alu_16_base_data as abd
from itertools import izip
import tables as tb
import os
filters = tb.Filters(complevel=5, complib='zlib')

'''
self.transcripts.dtype
Out[4]: dtype([('transcript_index', '<u4'), ('gene_region_index', '<u4'), ('gene_name', 'S10'),
 ('gene_symbol', 'S40'), ('strand', 'S1'), ('chrom', '<u2'), ('tx_start', '<u4'), ('tx_end', '<u4'),
 ('cds_start', '<u4'), ('cds_end', '<u4'), ('exon_count', '<u4')])

self.exons.dtype
Out[5]: dtype([('gene_region_index', '<u4'), ('transcript_index', '<u4'), ('chrom', '<u2'),
 ('exon_start', '<u4'), ('exon_end', '<u4'), ('strand', 'S1')])
'''

class transcript_data_cls(object) :
    
    data_folder = 'grch37_hg19_annotation_data'        
    file_name = 'all_known_gene_transcripts.h5'
    file_path = os.path.join(data_folder, file_name)
    transcript_table_name = 'all_known_gene_transcripts'
    exon_table_name = 'transcript_exons'
    gene_table_name = 'known_gene_regions'
    
    def __init__(self) :
        self.read_data()
        
    def read_data(self) :
        h5 = tb.open_file(self.file_path, 'r')
        transcript_table = getattr(h5.root, self.transcript_table_name)
        self.transcripts = transcript_table[:]
        exon_table = getattr(h5.root, self.exon_table_name)
        self.exons = exon_table[:]
        gene_table = getattr(h5.root, self.gene_table_name)
        self.genes = gene_table[:]
        h5.close()



class alu_strand_cg_base_cls(object) :        
    data_folder = 'grch37_hg19_alu_data'
    cg_count_dtype = np.dtype([('cg_count', np.uint32), ('count', np.uint32)])
    
    def __init__(self) :
        self.read_data()
        
    def read_data(self) :
        file_path = os.path.join(self.data_folder, self.file_name)
        h5 = tb.open_file(file_path, 'r')
        repeat_table = getattr(h5.root, self.alu_repeat_table_name)
        self.repeats = repeat_table[:]
        h5.close()

class pos_alu_cg_cls(alu_strand_cg_base_cls) :
    file_name = 'pos_strand_repeat_cg_data.h5'    
    alu_repeat_table_name = 'alu_plus_strand_repeats'
    
    def __init__(self) :
        alu_strand_cg_base_cls.__init__(self)
        
        
class neg_alu_cg_cls(alu_strand_cg_base_cls) :
    file_name = 'neg_strand_repeat_cg_data.h5'    
    alu_repeat_table_name = 'alu_minus_strand_repeats'
    
    def __init__(self) :
        alu_strand_cg_base_cls.__init__(self)


class cg_repeat_rdr_cls(object) :
    repeat_descr = [('index', '<u4'), ('chrom', '<u2'), ('start_pos', '<u4'), ('end_pos', '<u4'),
                    ('strand', '|S1'), ('repeat_name', '|S20'),  ('repeat_class', '|S20'), ('repeat_family', '|S20'),
                    ('sw_score', '<u2'), ('mismatches', '<u2'), ('deleted', '<u2'), ('inserted', '<u2'),
                    ('repeat_start', '<i2'), ('repeat_end', '<i2'), ('repeat_left', '<i2'), ('cg_count', '<u2'),
                    ('genome_region_index', np.uint32), ('genome_region_type', 'S1'), ('transcript_strand', 'S1')]    

    in_field_names = ('index', 'chrom', 'start_pos', 'end_pos', 'strand', 'repeat_name', 'repeat_class',
                      'repeat_family', 'sw_score', 'mismatches', 'deleted',  'inserted',  'repeat_start',
                      'repeat_end', 'repeat_left', 'cg_count')
    
    def build_out_repeats(self, in_repeats) :
        out_repeats = np.zeros(in_repeats.size, dtype=np.dtype(self.repeat_descr))
        for name in self.in_field_names :
            out_repeats[name] = in_repeats[name]
        return out_repeats
    
    def read_pos_repeats(self) :
        pro = pos_alu_cg_cls()
        pos_repeats = pro.repeats
        return self.build_out_repeats(pos_repeats)
    
    def read_neg_repeats(self) :
        nro = neg_alu_cg_cls()
        neg_repeats = nro.repeats
        return self.build_out_repeats(neg_repeats)
        

class genome_regions_cls(object) :
    chrom_region_dtype = np.dtype([('region_index', np.uint32), ('chrom', np.uint16), ('region_start', np.uint32),
                                   ('region_bound', np.uint32), ('region_type', 'S1'), ('transcript_strand', 'S1'),
                                   ('pos_alu', np.uint16), ('neg_alu', np.uint16)])
    genome_region_transcript_assoc_dtype = np.dtype([('genome_region_index', np.uint32), ('transcript_index', np.uint32)])
    chrom_region_none = (0, 0, 0, 0, '', '', 0, 0)
    def __init__(self) :
        self.read_data()
        
    def read_data(self) :
        trs = transcript_data_cls()
        self.transcripts = trs.transcripts
        self.exons = trs.exons
        self.repeat_rdr = cg_repeat_rdr_cls()
        self.genome_region_transcript_assoc = []

    def process_chrom_exons(self, chrom, chrom_exons) :
        chrom_regions = []
        exon_starts = chrom_exons['exon_start']
        exon_ends = chrom_exons['exon_end']
        exon_region_start = exon_starts[0]
        exon_region_bound = exon_ends[0]
        if exon_region_start > 1 :
            intron_region_start = 1
            intron_region_bound = exon_region_start
            region = (0,chrom, intron_region_start, intron_region_bound, 'o', '', 0, 0)
            chrom_regions.append(region)
        for exon_start, exon_end in izip(exon_starts[1:], exon_ends[1:]) :
            if exon_start > exon_region_bound :
                region = (0,chrom, exon_region_start, exon_region_bound, 'e', '', 0, 0)
                chrom_regions.append(region)
                intron_region_start = exon_region_bound
                intron_region_bound = exon_start
                region = (0,chrom, intron_region_start, intron_region_bound, 'o', '', 0, 0)
                chrom_regions.append(region)
                exon_region_start = exon_start
                exon_region_bound = exon_end
            else :
                if exon_end > exon_region_bound :
                    exon_region_bound = exon_end
        region = (0,chrom, exon_region_start, exon_region_bound, 'e', '', 0, 0)
        chrom_regions.append(region)
        region = (0,chrom, exon_region_bound, exon_region_bound, 'o', '', 0, 0)
        chrom_regions.append(region)
        chrom_regions = np.array(chrom_regions, dtype=self.chrom_region_dtype)
        return chrom_regions
                                                    
    def chrom_exons(self, chrom) :
            inds = self.exons['chrom'].searchsorted([chrom, chrom+1])
            exons = self.exons[inds[0]:inds[1]]
            return exons

    def process_exons(self) :
        self.exons.sort(order=['chrom', 'exon_start'])
        genome_regions = []
        region_none = np.array([self.chrom_region_none], dtype=self.chrom_region_dtype)
        genome_regions.append(region_none)
        for chrom in range(1,23) :
            chrom_exons = self.chrom_exons(chrom)
            chrom_regions = self.process_chrom_exons(chrom, chrom_exons)
            genome_regions.append(chrom_regions)
        self.genome_regions = np.concatenate(genome_regions)
        region_indexes = np.arange(self.genome_regions.size)
        self.genome_regions['region_index'] = region_indexes
    
        
    def update_transcript_regions(self, transcript, transcript_regions) :
        transcript_strand = transcript['strand']
        transcript_index = transcript['transcript_index']
        region_indexes = transcript_regions['region_index']
        for region_index in region_indexes :
            self.genome_region_transcript_assoc.append((region_index, transcript_index))
        m = transcript_regions['region_type'] == 'o'
        transcript_regions['region_type'][m] = 'i'
        region_strands = transcript_regions['transcript_strand']
        m = region_strands == ''
        region_strands[m] = transcript_strand
        m = region_strands != transcript_strand
        region_strands[m] = 'b'
        
    def process_chrom_transcripts(self, chrom, transcripts, regions) :
        transcript_starts = transcripts['tx_start']
        transcript_ends = transcripts['tx_end']
        region_starts = regions['region_start']
        inds_region_starts = region_starts.searchsorted(transcript_starts)
        inds_region_bounds = region_starts.searchsorted(transcript_ends)
        for i in range(transcripts.size) :
            transcript_regions = regions[inds_region_starts[i]:inds_region_bounds[i]]
            transcript = transcripts[i]
            self.update_transcript_regions(transcript, transcript_regions)
        

    def chrom_transcripts(self, chrom) :
        inds = self.transcripts['chrom'].searchsorted([chrom, chrom+1])
        transcripts = self.transcripts[inds[0]:inds[1]]
        return transcripts

    def chrom_regions(self, chrom) :
        inds = self.genome_regions['chrom'].searchsorted([chrom, chrom+1])
        regions = self.genome_regions[inds[0]:inds[1]]
        return regions
        
    def process_transcripts(self) :
        self.genome_regions.sort(order=['chrom', 'region_start', 'region_bound'])
        self.transcripts.sort(order=['chrom', 'tx_start', 'tx_end'])
        for chrom in range(1,23) :
            chrom_regions = self.chrom_regions(chrom)
            chrom_transcripts = self.chrom_transcripts(chrom)
            self.process_chrom_transcripts(chrom, chrom_transcripts, chrom_regions)
        self.genome_region_transcript_assoc = np.array(self.genome_region_transcript_assoc,
                                     dtype=self.genome_region_transcript_assoc_dtype)
        
        

    def process_chrom_repeats(self, chrom, chrom_regions, chrom_repeats, repeat_strand) :
        repeat_starts = chrom_repeats['start_pos']
        inds_repeat_region_start = repeat_starts.searchsorted(chrom_regions['region_start'])
        inds_repeat_region_bound = repeat_starts.searchsorted(chrom_regions['region_bound'])
        for i in range(chrom_regions.size) :
            region = chrom_regions[i]
            region_repeats = chrom_repeats[inds_repeat_region_start[i]:inds_repeat_region_bound[i]]
            region_repeat_count = region_repeats.size
            if region_repeat_count == 0 :
                continue
            region_repeats['genome_region_index'] = region['region_index']
            region_repeats['genome_region_type'] = region['region_type']
            region_repeats['transcript_strand'] = region['transcript_strand']
            if repeat_strand == '+' :
                region['pos_alu'] = region_repeat_count
            else :
                region['neg_alu'] = region_repeat_count
        assert i == chrom_regions.size - 1
        if region_repeat_count > 0 :
            last_region_bound = region['region_bound']
            last_region_repeat = region_repeats[-1]
            last_repeat_end_pos = last_region_repeat['end_pos']
            if last_repeat_end_pos > last_region_bound :
                region['region_bound'] = last_repeat_end_pos + 1
            

    def chrom_repeats(self, chrom, repeats) :
        inds = repeats['chrom'].searchsorted([chrom, chrom+1])
        return repeats[inds[0]:inds[1]]
        
    def process_repeats(self, repeats, repeat_strand) :
        for chrom in range(1, 23) :
            chrom_repeats = self.chrom_repeats(chrom, repeats)
            chrom_regions = self.chrom_regions(chrom)
            self.process_chrom_repeats(chrom, chrom_regions, chrom_repeats, repeat_strand)
        
    def process_pos_repeats(self) :
        self.pos_repeats = self.repeat_rdr.read_pos_repeats()
        self.process_repeats(self.pos_repeats, '+')

    def process_neg_repeats(self) :            
        self.neg_repeats = self.repeat_rdr.read_neg_repeats()
        self.process_repeats(self.neg_repeats, '-')
        
    def do_work(self) :
        self.process_exons()
        self.process_transcripts()
        self.process_pos_repeats()
        self.process_neg_repeats()


class region_data_writer_cls(object) :
    data_folder = 'grch37_hg19_alu_data'
    pos_repeat_file_name = 'pos_strand_repeat_cg_region_data.h5'    
    pos_alu_repeat_table_name = 'alu_plus_strand_repeats'
    neg_repeat_file_name = 'neg_strand_repeat_cg_region_data.h5'    
    neg_alu_repeat_table_name = 'alu_minus_strand_repeats'
    genome_region_file_name = 'genome_region_data.h5'
    genome_regions_table_name = 'genome_regions'
    genome_region_transcript_assoc_table_name = 'genome_region_transcript_associations'
    
    def __init__(self) :
        self.gro = genome_regions_cls()
        self.gro.do_work()
        
    def write_gro_pos_repeats(self) :
        file_path = os.path.join(self.data_folder, self.pos_repeat_file_name)
        h5_out = tb.open_file(file_path, 'w', filters=filters)
        pos_table = h5_out.create_table('/', self.pos_alu_repeat_table_name, description=self.gro.pos_repeats.dtype)
        pos_table.append(self.gro.pos_repeats)
        h5_out.close()
        
    def write_gro_neg_repeats(self) :
        file_path = os.path.join(self.data_folder, self.neg_repeat_file_name)
        h5_out = tb.open_file(file_path, 'w', filters=filters)
        neg_table = h5_out.create_table('/', self.neg_alu_repeat_table_name, description=self.gro.neg_repeats.dtype)
        neg_table.append(self.gro.neg_repeats)
        h5_out.close()
        
    def write_region_data(self) :
        file_name = self.genome_region_file_name
        file_path = os.path.join(self.data_folder, file_name)
        h5_out = tb.open_file(file_path, 'w', filters=filters)
        region_table = h5_out.create_table('/', self.genome_regions_table_name, description=self.gro.genome_regions.dtype)
        region_table.append(self.gro.genome_regions)
        assoc_table = h5_out.create_table('/', self.genome_region_transcript_assoc_table_name,
                                            description=self.gro.genome_region_transcript_assoc.dtype)
        assoc_table.append(self.gro.genome_region_transcript_assoc)
        h5_out.close()
        
    def do_work(self) :
        self.write_gro_pos_repeats()
        self.write_gro_neg_repeats()
        self.write_region_data()
        
'''        
rdwo = region_data_writer_cls()
rdwo.do_work()        
'''    
    














