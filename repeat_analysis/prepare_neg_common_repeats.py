# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
np.set_printoptions(precision=3, suppress=True)
import genomes_dnj_2.repeat_analysis_data.alu_16_base_data as abd
import genomes_dnj_2.repeat_analysis.offset_sequences as osq
import genomes_dnj_2.repeat_analysis.num_16_decode as n16d
import genomes_dnj_2.repeat_analysis.offset__cluster_matches as ocm
import genomes_dnj_2.repeat_analysis.num_partition_analysis as npa
import genomes_dnj_2.repeat_analysis.alu_sequence_alignments as align
import tables as tb
#filters = tb.Filters(complevel=5, complib='zlib')



num_16_26 = 3277099518
num_16_117 = 3338678528
num_16_245 =  3832882666
num_16_178 = 2323162274
num_16_282 = 1073741824


data_obj = abd.neg_alu_data_cls()
h5 = tb.open_file('neg_ris_most_common_num_16.h5', 'w')
ndo = data_obj.num_data_obj
osqo = osq.offset_sequence_cls(ndo)

a_20 = align.alignments[20]
onmo = ocm.offset_num_16_match_cls(a_20.base_offsets, osqo)
ris_neg_common_26 = onmo.ris_num_16(6, num_16_26)

h5.create_array('/', 'ris_neg_common_26', obj=ris_neg_common_26)


a_117 = align.alignments[117]
onmo = ocm.offset_num_16_match_cls(a_117.base_offsets, osqo)
ris_neg_common_117 = onmo.ris_num_16(0, num_16_117)


h5.create_array('/', 'ris_neg_common_117', obj= ris_neg_common_117)




a_178 = align.alignments[178]
onmo = ocm.offset_num_16_match_cls(a_178.base_offsets, osqo)
ris_neg_common_178 = onmo.ris_num_16(0, num_16_178)


h5.create_array('/', 'ris_neg_common_178', obj= ris_neg_common_178)



a_245 = align.alignments[245]
onmo = ocm.offset_num_16_match_cls(a_245.base_offsets, osqo)
ris_neg_common_245 = onmo.ris_num_16(0, num_16_245)


h5.create_array('/', 'ris_neg_common_245', obj= ris_neg_common_245)



a_282 = align.alignments[282]
onmo = ocm.offset_num_16_match_cls(a_282.base_offsets, osqo)
ris_neg_common_282 = onmo.ris_num_16(0, num_16_282)


h5.create_array('/', 'ris_neg_common_282', obj= ris_neg_common_282)


h5.close()


'''
ris_pos_common_26.size
Out[2]: 189341

ris_pos_common_117.size
Out[3]: 166724

ris_pos_common_245.size
Out[4]: 157218

ris_pos_common_178.size
Out[5]: 150835

ris_pos_common_282.size
Out[6]: 101327

'''







