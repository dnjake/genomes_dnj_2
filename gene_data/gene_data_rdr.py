# -*- coding: utf-8 -*-

import tables as tb
import os

s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)

file_name = 'chrom_filtered_genes.h5'
file_path = os.path.join(mod_dir, file_name)

def table_name(chrom) :
    return 'chrom_' + str(chrom) + '_genes'


class genes_rdr_cls(object) :
    def __init__(self, chrom) :
        h5 = tb.open_file(file_path, 'r')
        table = getattr(h5.root, table_name(chrom))
        self.genes = table[:]
        h5.close()
        
    def genes_in_interval(self, first, last) :
        ind_first = self.genes['tx_end'].searchsorted(first)
        ind_last = self.genes['tx_start'].searchsorted(last)
        return self.genes[ind_first:ind_last]