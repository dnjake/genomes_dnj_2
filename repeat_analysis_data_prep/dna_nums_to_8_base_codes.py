# -*- coding: utf-8 -*-


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')


class build_dna_seq_nums_by_pos_cls(object) :
    '''
    This class aggregates 8 single base numerical
    codes into one 16 bit number.  A 16 bit number
    is created for each position in the input data
    that is followed by at least 7 more positions.
    
    A valid mask is created that is true for all
    positions where all of the 8 input numbers represent
    valid DNA bases.
    
    Invalid position have seq_num 0.
    But seq_num 0 is also one of the most common valid values
    
    Input values are 0, 1, 2, or 3 for valid DNA bases
    Invalid input positions have value 128
    '''
    def __init__(self, by_pos_single_base_data) :
        self.by_pos_single_base_data = by_pos_single_base_data
        
    def build_by_pos_sequence_numbers(self) :
        cd = self.by_pos_single_base_data
        cds = []
        for i in range(7) :
            cds.append(cd[i:-(7-i)])
        cds.append(cd[7:])
        cds_size = cds[0].size        
        seq_num_data = np.zeros(cds_size, dtype=np.uint16)
        seq_num_mask = np.zeros(cds_size, dtype=np.bool)                
        seq_num_mask = cds[0] == 128
        for s in cds[1:] :
            seq_num_mask = np.logical_or(seq_num_mask, s==128)        
        seq_num_mask = np.logical_not(seq_num_mask)
        self.by_pos_valid_mask = seq_num_mask        
        s = cds[0]        
        seq_num_data[seq_num_mask] = s[seq_num_mask]
        for s in cds[1:] :
            seq_num_data[seq_num_mask] = np.left_shift(seq_num_data[seq_num_mask], 2)
            seq_num_data[seq_num_mask] = np.bitwise_or(seq_num_data[seq_num_mask], s[seq_num_mask])
        self.by_pos_sequence_numbers = seq_num_data
    


class chrom_base_to_sequence_nums_cls(object) :
    base_numbers_folder = 'grch37_hg19_chrom_dna_numbers'
    base_numbers_file_name_end = '_base_numbers_by_pos.h5'
    base_numbers_carray_name_end = '_base_numbers_by_pos'
    sequence_numbers_folder = 'grch37_hg19_chrom_dna_8_base_codes'
    sequence_numbers_file_name_end = '_8_base_sequence_numbers_by_pos.h5'
    sequence_numbers_carray_name_end = '_8_base_sequence_numbers_by_pos'
    valid_carray_name_end = '_valid_pos_mask'
    def __init__(self, chrom) :
        self.chrom = chrom
        self.chrom_start_str = 'chr_' + str(self.chrom)
    
    def read_base_numbers(self) :
        file_path = self.base_numbers_folder + '/' + self.chrom_start_str + self.base_numbers_file_name_end
        carray_name = self.chrom_start_str + self.base_numbers_carray_name_end
        h5_in = tb.open_file(file_path, 'r')
        carray = getattr(h5_in.root, carray_name)
        self.dna_base_numbers = carray[:]
        h5_in.close()
        
    def sequences_from_bases(self) :
        bsno = build_dna_seq_nums_by_pos_cls(self.dna_base_numbers)
        bsno.build_by_pos_sequence_numbers()
        self.dna_sequence_numbers = bsno.by_pos_sequence_numbers
        self.valid_mask = bsno.by_pos_valid_mask
        
    def write_sequence_numbers(self) :
        file_path = self.sequence_numbers_folder + '/' + self.chrom_start_str + self.sequence_numbers_file_name_end
        seq_num_carray_name = self.chrom_start_str + self.sequence_numbers_carray_name_end
        valid_mask_name = self.chrom_start_str + self.valid_carray_name_end
        h5_out = tb.open_file(file_path, 'w', filters=filters)
        h5_out.create_carray('/', seq_num_carray_name, obj=self.dna_sequence_numbers)
        h5_out.create_carray('/', valid_mask_name, obj=self.valid_mask)
        h5_out.close()
        
    def do_transform(self) :
        self.read_base_numbers()
        self.sequences_from_bases()
        self.write_sequence_numbers()

'''
chroms = range(1, 23)

for chrom in chroms :
    print ('transforming chrom ' + str(chrom))
    snfbo = chrom_base_to_sequence_nums_cls(chrom)
    snfbo.do_transform()
    
print('done')
'''



