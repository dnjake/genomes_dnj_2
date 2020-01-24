# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')


class updater_cls(object) :


    def update_pos_alu(self) :
        in_alu_path = os.path.join('grch37_hg19_alu_data', 'pos_alu_dna_num_16_initial.h5')
        h5_in = tb.open_file(in_alu_path, 'r')
        in_alu_table = getattr(h5_in.root, 'pos_alu_dna_num_16_data')
        pos_alu_data = in_alu_table[:]
        in_repeat_path = os.path.join('grch37_hg19_alu_data', 'pos_strand_repeat_cg_region_data.h5')
        h5_repeats = tb.open_file(in_repeat_path, 'r')
        in_repeat_table = getattr(h5_repeats.root, 'alu_plus_strand_repeats')
        in_repeat_data = in_repeat_table[:]
        out_path = 'pos_alu_dna_num_16.h5'
        h5_out = tb.open_file(out_path, 'w', filters=filters)
        out_alu_table = h5_out.create_table('/', 'pos_alu_dna_num_16_data', description=pos_alu_data.dtype)
        out_alu_table.append(pos_alu_data)
        out_repeat_table = h5_out.create_table('/', 'alu_plus_strand_repeats', description=in_repeat_data.dtype)
        out_repeat_table.append(in_repeat_data)
        h5_out.close()
        h5_in.close()
        h5_repeats.close()


    def update_neg_alu(self) :
        in_alu_path = os.path.join('grch37_hg19_alu_data', 'neg_alu_dna_num_16_initial.h5')
        h5_in = tb.open_file(in_alu_path, 'r')
        in_alu_table = getattr(h5_in.root, 'neg_alu_dna_num_16_data')
        neg_alu_data = in_alu_table[:]
        in_repeat_path = os.path.join('grch37_hg19_alu_data', 'neg_strand_repeat_cg_region_data.h5')
        h5_repeats = tb.open_file(in_repeat_path, 'r')
        in_repeat_table = getattr(h5_repeats.root, 'alu_minus_strand_repeats')
        in_repeat_data = in_repeat_table[:]
        out_path = 'neg_alu_dna_num_16.h5'
        h5_out = tb.open_file(out_path, 'w', filters=filters)
        out_alu_table = h5_out.create_table('/', 'neg_alu_dna_num_16_data', description=neg_alu_data.dtype)
        out_alu_table.append(neg_alu_data)
        out_repeat_table = h5_out.create_table('/', 'alu_minus_strand_repeats', description=in_repeat_data.dtype)
        out_repeat_table.append(in_repeat_data)
        h5_out.close()
        h5_in.close()
        h5_repeats.close()



'''
updo = updater_cls()
updo.update_pos_alu()
updo.update_neg_alu()
'''

















