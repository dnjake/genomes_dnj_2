# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')



'''
This just writes numbers for each letter
and writes 128 in cases where the letter is invalid
'''


class chrom_dna_rdr_cls(object) :
    letters_folder = 'grch37_hg19_chrom_dna_letters'
    file_name_start = 'grch37_chr'
    file_name_end = '_fa.h5'
    table_prefix = 'grch37_chr_'
    def __init__(self, chrom) :
        chrom_str = str(chrom)
        self.file_name = self.file_name_start + chrom_str + self.file_name_end
        self.file_path = os.path.join(self.letters_folder, self.file_name)
        self.table_name = self.table_prefix + chrom_str
        self.read_chrom_dna()
        
    def read_chrom_dna(self) :
        h5_in = tb.open_file(self.file_path, 'r')
        data_table = getattr(h5_in.root, self.table_name)
        self.dna_data = data_table[:]
        h5_in.close()

class dna_sequence_letters_to_numbers_cls(object) :
    dna_base_codes = (('A', 0), ('C', 1), ('G', 2), ('T', 3))
    def __init__(self, dna_letters) :
        self.dna_letters = dna_letters
        
    def do_transform(self) :
        dna_letters = self.dna_letters
        base_data = np.empty(dna_letters.size, dtype=np.uint8)
        base_data[:] = 128
        for letter, number in self.dna_base_codes :
            m = dna_letters == letter
            base_data[m] = number
        self.dna_numbers = base_data

class chrom_letters_to_numbers_transform_cls(object) :
    chrom_numbers_folder = 'grch37_hg19_chrom_dna_numbers'
    numbers_file_name_end = '_base_numbers_by_pos.h5'
    carray_name_end = '_base_numbers_by_pos'
    def __init__(self, chrom) :
        self.chrom = chrom
        
    def read_chrom_letters(self) :
        fro = chrom_dna_rdr_cls(self.chrom)
        fro.read_chrom_dna()
        self.chrom_dna_letters = fro.dna_data
        
    def transform_letters(self):
        tlo = dna_sequence_letters_to_numbers_cls(self.chrom_dna_letters)
        tlo.do_transform()
        self.chrom_dna_numbers = tlo.dna_numbers
        
    def write_dna_numbers(self) :
        chrom = str(self.chrom)
        chrom_start = 'chr_' + chrom
        file_path = self.chrom_numbers_folder + '/' + chrom_start + self.numbers_file_name_end
        carray_name = chrom_start + self.carray_name_end
        h5 = tb.open_file(file_path, 'w', filters=filters)
        h5.create_carray('/', carray_name, obj=self.chrom_dna_numbers)
        h5.close()
        
    def do_chrom_transform(self) :
        self.read_chrom_letters()
        self.transform_letters()
        self.write_dna_numbers()


'''
chroms = range(1, 23)

for chrom in chroms :
    print ('transforming chrom ' + str(chrom))
    lnto = chrom_letters_to_numbers_transform_cls(chrom)
    lnto.do_chrom_transform()
    
print('done')
'''
