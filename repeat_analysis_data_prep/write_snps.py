# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import numpy as np
import tables as tb

import genomes_dnj_2.repeat_analysis_data.alu_16_base_data as abd
import genomes_dnj_2.autosome_snp_data.chrom_snp_data_rdr as csd
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d
filters = tb.Filters(complevel=5, complib='zlib')

class data_writer_cls(object) :
    alu_sequence_data_folder = 'grch37_hg19_alu_data'
    repeat_snp_dtype = np.dtype([('chrom', np.uint16), ('alu_offset', np.uint16), ('snp_index', '<u4'), ('pos', '<u4'),
                                 ('repeat_index', np.uint32) ])
    snp_data_dtype = np.dtype([('snp_index', '<u4'), ('pos', '<u4')])
    snp_stats_dtype = np.dtype([('chrom', np.uint16), ('all_snps', np.uint32), ('alu_snps', np.uint32)])

    def __init__(self, repeats) :
        self.repeats = repeats
    
    def read_chrom_snps(self, chrom) :
        rdr = csd.chrom_snp_data_tables_cls(chrom)
        snps = rdr.snp_data_table[:]
        rdr.close()
        snp_data = np.zeros(snps.size, dtype=self.snp_data_dtype)
        snp_data['snp_index'] = snps['snp_index']
        snp_data['pos'] = snps['pos']
        return snp_data
    
    def read_chrom_repeats(self, chrom) :
        bound = chrom + 1
        inds = self.repeats['chrom'].searchsorted([chrom, bound])
        return self.repeats[inds[0]:inds[1]]
    
    def chrom_alu_snps_repeat_indexes(self, snps, repeats) :    
        irs = repeats['end_pos'].searchsorted(snps['pos'])
        m = irs < repeats.size
        irs = irs[m]
        snps = snps[m]
        m = snps['pos'] > repeats['start_pos'][irs]
        irs = irs[m]
        snps = snps[m]
        return snps, irs
    
    def chrom_alu_snp_data(self, chrom) :
        chrom_snps = self.read_chrom_snps(chrom)
        chrom_repeats = self.read_chrom_repeats(chrom)
        chrom_alu_snps, inds_chrom_repeats = self.chrom_alu_snps_repeat_indexes(chrom_snps, chrom_repeats)
        snps = chrom_alu_snps
        sris = inds_chrom_repeats
        oda = np.zeros(snps.size, dtype=self.repeat_snp_dtype)
        oda['chrom'] = chrom
        oda['snp_index'] = snps['snp_index']
        oda['pos'] = snps['pos']
        oda['repeat_index'] = chrom_repeats['index'][sris]
        self.snp_alu_offsets(oda, chrom_repeats, sris)
        stats_data = chrom, chrom_snps.size, chrom_alu_snps.size
        return oda, stats_data
    
    def process_chroms(self) :
        out_data = []
        out_stats = []
        for chrom in range(1, 23) :
            print('processing', chrom)
            chrom_snp_data, chrom_stats_data = self.chrom_alu_snp_data(chrom)
            out_data.append(chrom_snp_data)
            out_stats.append(chrom_stats_data)
        self.snp_data = np.concatenate(out_data)
        self.alu_snp_stats = np.array(out_stats, dtype=self.snp_stats_dtype)
        
    def write_data(self) :
        file_path = os.path.join(self.alu_sequence_data_folder, self.file_name)
        h5 = tb.open_file(file_path, 'w', filters=filters)
        snp_data_table = h5.create_table('/', self.snp_data_table_name, description=self.snp_data.dtype)
        snp_data_table.append(self.snp_data)
        repeat_table = h5.create_table('/', self.alu_repeat_table_name, description=self.repeats.dtype)
        repeat_table.append(self.repeats)
        h5.close()
        
    def do_work(self) :
        self.process_chroms()
        self.write_data()

class pos_strand_data_writer_cls(data_writer_cls) :
    file_name = 'pos_alu_repeat_snps.h5'
    snp_data_table_name = 'pos_alu_repeat_snps'
    alu_repeat_table_name = 'alu_plus_strand_repeats'
    
    def __init__(self) :
        pos_repeats_obj = abd.pos_alu_cg_cls()
        pos_repeats = pos_repeats_obj.repeats
        data_writer_cls.__init__(self, pos_repeats)
        
    '''
    The problem is that the start position given in the repeat data is 1 less
    than the start position given on the browser.
    The reason is that the dna positions for the browser and my snp positions
    are 1 based.
    
    So the actual start position is 1 more than the data start position
    and the alu offset at position is at least 1
    '''

    def snp_alu_offsets(self, snp_data, repeats, snp_repeat_indexes) :
        sris = snp_repeat_indexes        
        alu_start_offsets = repeats['repeat_start'][sris]
        alu_start_pos = repeats['start_pos'][sris] + 1
        snp_pos_offsets = snp_data['pos'] - alu_start_pos
        snp_data['alu_offset'] = alu_start_offsets + snp_pos_offsets
        

class neg_strand_data_writer_cls(data_writer_cls) :
    file_name = 'neg_alu_repeat_snps.h5'
    snp_data_table_name = 'neg_alu_repeat_snps'
    alu_repeat_table_name = 'alu_minus_strand_repeats'
    
    def __init__(self) :
        neg_repeats_obj = abd.neg_alu_cg_cls()
        neg_repeats = neg_repeats_obj.repeats
        data_writer_cls.__init__(self, neg_repeats)
    '''
        data_size = alu_first_seqs.size - 15
        alu_data = np.zeros(data_size, self.alu_seq_dtype)
        alu_data['alu_offset'] = np.arange(data_size, dtype=np.uint16)
        alu_data['alu_offset'] = data_size - alu_data['alu_offset'] -1
        alu_data = np.zeros(data_size, self.alu_seq_dtype)
        alu_data['chrom'] = chrom
        alu_data['alu_offset'] = np.arange(data_size, dtype=np.uint16)
        alu_data['alu_offset'] = data_size - alu_data['alu_offset'] -1
        alu_data['alu_offset'] += alu_start_offset

    '''        

    def snp_alu_offsets(self, snp_data, repeats, snp_repeat_indexes) :
        sris = snp_repeat_indexes        
        alu_bounds = repeats['end_pos'][sris]
        snp_pos_offsets = alu_bounds - snp_data['pos']
        snp_data['alu_offset'] = snp_pos_offsets + repeats['repeat_left'][sris]




'''   
dwo = neg_strand_data_writer_cls()
dwo.do_work()
'''

'''
possibly there are no snps in the first offset
dwo.repeats[138]
Out[34]: (1311, 1, 830669, 830888, '-', 'aluyk4', 'sine', 'alu', 1848, 59, 5, 0, -92, 220, 1, 12)

       (1, 150, 1550, 830739, 1311), (1, 111, 1551, 830778, 1311),
       (1, 107, 1552, 830782, 1311), (1, 106, 1553, 830783, 1311),
       (1, 104, 1554, 830785, 1311), (1,  52, 1555, 830837, 1311),
       (1,   8, 1556, 830881, 1311), (1, 278, 1587, 834430, 1318),
       
830888 - 830739       

'''

'''
dwo.alu_snp_stats
Out[2]: 
array([( 1, 1443793, 108498), ( 2, 1547546,  92197),
       ( 3, 1301411,  75232), ( 4, 1329202,  69308),
       ( 5, 1170056,  64086), ( 6, 1185982,  70856),
       ( 7, 1083764,  76706), ( 8, 1026789,  58757),
       ( 9,  800795,  54442), (10,  919273,  63552),
       (11,  906615,  52911), (12,  877190,  64592),
       (13,  655907,  35987), (14,  597788,  41446),
       (15,  541959,  41294), (16,  587078,  54895),
       (17,  504449,  59252), (18,  518146,  30040),
       (19,  422925,  63606), (20,  403986,  30925),
       (21,  255365,  15747), (22,  256408,  26960)],
      dtype=[('chrom', '<u2'), ('all_snps', '<u4'), ('alu_snps', '<u4')])

d = dwo.alu_snp_stats

snp_count = d['all_snps'].sum()

alu_count = d['alu_snps'].sum()

float(alu_count)/float(snp_count)
Out[6]: 0.0682406119796403

snp_count
Out[7]: 18336427

alu_count
Out[8]: 1251289


51880 is the actual position of the last nucleotide in the repeat
it is the position that should have an alu offset of 1 in this case
self.repeats[5]
Out[5]: (64, 1, 51584, 51880, '+', 'aluy', 'sine', 'alu', 2368, 74, 3, 0, 1, 297, -14, 29)

array([(1, 120,   34,  51762,   64), (1, 117,   35,  51765,   64),
       (1,  70,  107,  81032,  110), (1,  25,  151,  91190,  128),
       (1,  41,  206, 237505,  344), (1, 233,  208, 247792,  364),
       (1, 247,  226, 405010,  585), (1, 246,  227, 405011,  585),
       (1, 305,  237, 526840,  745), (1, 256,  238, 526889,  745),
       (1, 184,  429, 583483,  857), (1,  74,  463, 666249, 1012),
       (1, 141,  471, 669590, 1021), (1, 190,  480, 672940, 1030),
       (1, 169,  481, 672961, 1030), (1,  37,  489, 684747, 1070),
       (1,   9,  497, 701131, 1114), (1,  57,  507, 705452, 1121),
       (1, 113,  524, 710195, 1128), (1, 280,  525, 710332, 1129),
       (1, 222,  526, 710390, 1129), (1, 194,  527, 710418, 1129),
       (1,  76,  528, 710536, 1129), (1,  37,  529, 710575, 1129),
       (1,  98,  539, 713337, 1137), (1, 180,  548, 714832, 1141),
       (1, 279,  549, 715074, 1142), (1, 213,  550, 715140, 1142),
       (1, 211,  551, 715142, 1142), (1, 197,  553, 715241, 1143),
       (1, 173,  554, 715265, 1143), (1, 167,  555, 715271, 1143),
       (1,  71,  556, 715367, 1143), (1, 137,  599, 723608, 1167),
       (1,   3,  600, 723742, 1167), (1, 259,  720, 734162, 1187),
       (1, 154,  721, 734267, 1187), (1,  89,  722, 734332, 1187),
       (1,  83,  723, 734338, 1187), (1,  81,  724, 734340, 1187),
       (1,  72,  725, 734349, 1187), (1,  38,  726, 734383, 1187),
       (1, 289,  757, 738539, 1194), (1,  72,  758, 738756, 1194),
       (1, 228, 1525, 827610, 1301), (1, 171, 1531, 828468, 1303),
       (1, 100, 1532, 828539, 1303), (1, 301, 1535, 828883, 1305),
       (1, 166, 1536, 829018, 1305), (1,  78, 1537, 829106, 1305)],
      dtype=[('chrom', '<u2'), ('alu_offset', '<u2'), ('snp_index', '<u4'), ('pos', '<u4'), ('repeat_index', '<u4')])

snps.dtype
Out[7]: dtype([('snp_index', '<u4'), ('chrom', '<u2'), ('ref', 'S1'), ('alt', 'S1'), ('id', 'S11'), 
('not_expressed_is_variant', 'u1'), ('pos', '<u4'), ('all_count', '<u2'), ('afr_obs_to_pred', '<f4'),
 ('afx_obs_to_pred', '<f4'), ('amr_obs_to_pred', '<f4'), ('eas_obs_to_pred', '<f4'), ('eur_obs_to_pred', '<f4'),
 ('sas_obs_to_pred', '<f4'), ('sax_obs_to_pred', '<f4'), ('acb_count', '<u2'), ('asw_count', '<u2'),
 ('beb_count', '<u2'), ('cdx_count', '<u2'), ('ceu_count', '<u2'), ('chb_count', '<u2'), ('chs_count', '<u2'),
 ('clm_count', '<u2'), ('esn_count', '<u2'), ('fin_count', '<u2'), ('gbr_count', '<u2'), ('gih_count', '<u2'),
 ('gwd_count', '<u2'), ('ibs_count', '<u2'), ('itu_count', '<u2'), ('jpt_count', '<u2'), ('khv_count', '<u2'),
 ('lwk_count', '<u2'), ('msl_count', '<u2'), ('mxl_count', '<u2'), ('pel_count', '<u2'), ('pjl_count', '<u2'),
 ('pur_count', '<u2'), ('stu_count', '<u2'), ('tsi_count', '<u2'), ('yri_count', '<u2')])
    
    
 dtype=[('index', '<u4'), ('chrom', '<u2'), ('start_pos', '<u4'), ('end_pos', '<u4'), ('strand', 'S1'),
 ('repeat_name', 'S20'), ('repeat_class', 'S20'), ('repeat_family', 'S20'), ('sw_score', '<u2'),
 ('mismatches', '<u2'), ('deleted', '<u2'), ('inserted', '<u2'), ('repeat_start', '<i2'),
 ('repeat_end', '<i2'), ('repeat_left', '<i2'), ('cg_count', '<u2')])    

'''