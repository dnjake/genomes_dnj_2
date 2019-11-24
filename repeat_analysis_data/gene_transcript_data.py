# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import tables as tb


s = __file__
mod_path = os.path.abspath(s)
mod_dir = os.path.dirname(mod_path)


class transcript_data_cls(object) :
    
    data_folder = 'grch37_hg19_annotation_data'        
    file_name = 'all_known_gene_transcripts.h5'
    local_file_path = os.path.join(data_folder, file_name)
    file_path = os.path.join(mod_dir, local_file_path)
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
   