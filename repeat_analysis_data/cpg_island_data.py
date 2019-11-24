# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import tables as tb


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)


class all_cpg_island_cls(object) :
    data_folder = 'grch37_hg19_annotation_data'
    file_name = 'cpg_all_islands.h5'
    local_file_path = os.path.join(data_folder, file_name)
    file_path = os.path.join(mod_dir, local_file_path)
    table_name = 'cpg_all_islands'
    
    def __init__(self) :
        self.read_data()

    def read_data(self) :
        h5 = tb.open_file(self.file_path, 'r')
        data_table = getattr(h5.root, self.table_name)
        self.all_cpg_islands = data_table[:]
        h5.close()
    
        
class masked_cpg_island_cls(object) :
    data_folder = 'grch37_hg19_annotation_data'
    file_name = 'cpg_masked_islands.h5'
    local_file_path = os.path.join(data_folder, file_name)
    file_path = os.path.join(mod_dir, local_file_path)
    table_name = 'cpg_masked_islands'

    def __init__(self) :
        self.read_data()

    def read_data(self) :
        h5 = tb.open_file(self.file_path, 'r')
        data_table = getattr(h5.root, self.table_name)
        self.masked_cpg_islands = data_table[:]
        h5.close()









    