# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

class chrom_dna_cls(object) :
    letters_folder = 'grch37_hg19_chrom_dna_letters'
    file_name_start = 'grch37_chr'
    file_name_end = '_fa.h5'
    array_name_start = 'grch37_chr_'

    def __init__(self, chrom) :
        self.chrom = chrom
        self.chrom_str = str(chrom)
        
    def read_file(self) :
        file_name = self.file_name_start + self.chrom_str + self.file_name_end
        local_path = os.path.join(self.letters_folder, file_name)
        file_path = os.path.join(mod_dir, local_path)
        array_name = self.array_name_start + self.chrom_str
        h5 = tb.open_file(file_path, 'r')
        array = getattr(h5.root, array_name)
        self.dna = array[:]
        h5.close()
        
        

        