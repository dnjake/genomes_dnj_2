# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')
import gzip

'''
source for grch37 chromosome data files is 
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/

This directory contains the Feb. 2009 assembly of the human genome (hg19,
GRCh37 Genome Reference Consortium Human Reference 37 (GCA_000001405.1))
in one gzip-compressed FASTA file per chromosome.
'''



class references_to_h5_cls(object) :
    reference_folder = 'grch37_hg19_ncbi_reference'
    in_file_name_start = 'chr'
    in_file_name_end = '.fa.gz'
    letters_folder = 'grch37_hg19_chrom_dna_letters'
    out_file_name_start = 'grch37_chr'
    out_file_name_end = '_fa.h5'
    out_table_name_start = 'grch37_chr_'    
    S1 = tb.StringAtom(itemsize=1)

    def convert_chrom(self, chrom) :
        in_file_name = self.in_file_name_start + str(chrom) + self.in_file_name_end
        in_file_path = os.path.join(self.reference_folder, in_file_name)
        statv = os.stat(in_file_path)
        erows = statv.st_size        
        f = gzip.open(in_file_path, 'r')
        data_lines = f.readlines()
        f.close()
        out_file_name = self.out_file_name_start + str(chrom) + self.out_file_name_end
        out_file_path = os.path.join(self.letters_folder, out_file_name)
        out_table_name = self.out_table_name_start + str(chrom)
        h5 = tb.open_file(out_file_path, 'w', filters=filters)
        data_array = h5.create_earray('/', out_table_name, self.S1, (0,), expectedrows=erows)
        data_array.append([''])
        for l in data_lines :
            if l[0] == '>' :
                continue            
            s = l.rstrip()
            sa = np.frombuffer(s,dtype='S1')
            data_array.append(sa)
        h5.close()
        
    def do_chroms(self) :
        for chrom in range(1, 23) :
            print ('transforming chrom ' + str(chrom))
            self.convert_chrom(chrom)
                
'''        
rtho = references_to_h5_cls()
rtho.do_chroms()
'''