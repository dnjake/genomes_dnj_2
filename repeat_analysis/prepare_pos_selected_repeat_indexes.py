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


data_obj = abd.pos_alu_data_cls()
h5 = tb.open_file('pos_selected_ris.h5', 'w')
ndo = data_obj.num_data_obj
osqo = osq.offset_sequence_cls(ndo)

aro_w = ocm.alignment_run_cls(align.alignments[48], osqo, ndo)
align_offset = aro_w.alu_offset

ocmo = aro_w.ocmo
bmos = ocmo.base_offset_matches
cmso = ocm.cluster_seq_match_ris_cls(bmos, ndo)
align_68 = align.alignments[68]
patterns_68 = align_68.patterns

with_62_63_offset_seqs = []
offset_delta_62_63 = 68 - align_offset

no_62_63_offset_seqs = []
offset_delta_no_62_63 = offset_delta_62_63 - 2

for p in patterns_68 :
    offset = p.delta_offset + offset_delta_62_63
    with_62_63_offset_seqs.append((offset, p.match_seqs))
    offset -= 2
    no_62_63_offset_seqs.append((offset, p.match_seqs))
    
ris_with_62_63 = cmso.or_pattern_match_ris(with_62_63_offset_seqs)
ris_no_62_63 = cmso.or_pattern_match_ris(no_62_63_offset_seqs)
matched_ris = np.union1d(ris_with_62_63, ris_no_62_63)
ris_aligned_48_no_match_68 = np.setdiff1d(aro_w.aligned_ris, matched_ris)

h5.create_array('/', 'ris_standard_62_63', obj=ris_with_62_63)
h5.create_array('/', 'ris_deleted_62_63', obj=ris_no_62_63)
h5.create_array('/', 'ris_aligned_48_no_match_68', obj=ris_aligned_48_no_match_68)

# also want ag at offset 57
o_57 = 57 - align_offset
ag_offsqs = [(o_57, ['AG'])]
ris_ag_57 = cmso.and_pattern_match_ris(ag_offsqs)
ris_aligned_48_not_ag_57 = np.setdiff1d(aro_w.aligned_ris, ris_ag_57)
h5.create_array('/', 'ris_ag_57', obj=ris_ag_57)
h5.create_array('/', 'ris_aligned_48_not_ag_57', obj=ris_aligned_48_not_ag_57)

# also want ta at offset 100

aro_w = ocm.alignment_run_cls(align.alignments[80], osqo, ndo)
align_offset = aro_w.alu_offset

ocmo = aro_w.ocmo
bmos = ocmo.base_offset_matches
cmso = ocm.cluster_seq_match_ris_cls(bmos, ndo)

o_100 = 100 - align_offset
ta_offsqs = [(o_100, ['TA'])]
ris_ta_100 = cmso.and_pattern_match_ris(ta_offsqs)
ris_aligned_80_not_ta_100 = np.setdiff1d(aro_w.aligned_ris, ris_ta_100)
h5.create_array('/', 'ris_ta_100', obj=ris_ta_100)
h5.create_array('/', 'ris_aligned_80_not_ta_100', obj=ris_aligned_80_not_ta_100)


# taat tagt at offset 161

aro_w = ocm.alignment_run_cls(align.alignments[159], osqo, ndo)
align_offset = aro_w.alu_offset

ocmo = aro_w.ocmo
bmos = ocmo.base_offset_matches
cmso = ocm.cluster_seq_match_ris_cls(bmos, ndo)

o_161 = 161 - align_offset
taat_offsqs = [(o_161, ['TAAT'])]
ris_taat_161 = cmso.and_pattern_match_ris(taat_offsqs)
tagt_offsqs = [(o_161, ['TAGT'])]
ris_tagt_161 = cmso.and_pattern_match_ris(tagt_offsqs)
matched_ris = np.union1d(ris_taat_161, ris_tagt_161)
ris_aligned_159_no_match_161 = np.setdiff1d(aro_w.aligned_ris, matched_ris)
h5.create_array('/', 'ris_taat_161', obj=ris_taat_161)
h5.create_array('/', 'ris_tagt_161', obj=ris_tagt_161)
h5.create_array('/', 'ris_aligned_159_no_match_161', obj=ris_aligned_159_no_match_161)

aro_w = ocm.alignment_run_cls(align.alignments[216], osqo, ndo)
align_offset = aro_w.alu_offset

ocmo = aro_w.ocmo
bmos = ocmo.base_offset_matches
cmso = ocm.cluster_seq_match_ris_cls(bmos, ndo)



o_197 = 197 - align_offset
o_199 = o_197 + 2

p_gg = [(o_197, ['GG']), (o_199, ['CG', 'CA', 'TG'])]
p_ct = [(o_199, ['CT', 'CC', 'TT'])]

ris_gg = cmso.and_pattern_match_ris(p_gg)
ris_ct = cmso.and_pattern_match_ris(p_ct)
matched_ris = np.union1d(ris_gg, ris_ct)
not_matched_ris = np.setdiff1d(aro_w.aligned_ris, matched_ris)


h5.create_array('/', 'ris_197_gg', obj=ris_gg)
h5.create_array('/', 'ris_199_ct', obj=ris_ct)
h5.create_array('/', 'ris_aligned_216_no_match_197', obj=not_matched_ris)


aro_w = ocm.alignment_run_cls(align.alignments[139], osqo, ndo)
align_offset = aro_w.alu_offset

ocmo = aro_w.ocmo
bmos = ocmo.base_offset_matches
cmso = ocm.cluster_seq_match_ris_cls(bmos, ndo)

o_153 = 153 - align_offset

p_gg = [(o_153, ['GG'])]

ris_gg = cmso.and_pattern_match_ris(p_gg)
matched_ris = ris_gg
not_matched_ris = np.setdiff1d(aro_w.aligned_ris, matched_ris)

h5.create_array('/', 'ris_153_gg', obj=ris_gg)
h5.create_array('/', 'ris_aligned_139_no_match_153', obj=not_matched_ris)


aro_w = ocm.alignment_run_cls(align.alignments[231], osqo, ndo)
align_offset = aro_w.alu_offset
ocmo = aro_w.ocmo
bmos = ocmo.base_offset_matches
cmso = ocm.cluster_seq_match_ris_cls(bmos, ndo)

seqs_30 = ['TGACA', 'CAACA', 'CGACA', 'AGACA', 'TAACA', 'CAATA', 'GGACA', 'TGATA', 'CCACA',
           'CTACA', 'TGGCA', 'CAGCA', 'CGATA', 'TGAAA', 'CAAAA', 'AAACA', 'CGAAA', 'CAAGA', 'TGATG', 
           'TGAGA', 'TGACG', 'CGGCA', 'CAACG', 'TGACT', 'GAACA',  'CAATG', 'CGACG', 'CGAGA']


seqs_35 = ['GAGTG', 'GAGCA', 'GAGCG', 'GAGTA', 'GAGGG', 'GAGAG', 'GAATG', 'GAGCC', 'AAGTG',
           'AAGCA', 'GAACA', 'CAGCA', 'CAGTG', 'AAGCG', 'CAGCG', 'GAACG', 'GAGCT', 'GGGTG', 
           'GGGCA', 'GGGCG', 'GAGAA', 'GTGTG', 'GAGGA', 'GACTG', 'GTGCG', 'GTGCA', 'GACCA']


offsqs_standard = [(30, seqs_30), (35, seqs_35)]
ris_standard_257 = cmso.and_pattern_match_ris(offsqs_standard)

h5.create_array('/', 'ris_standard_257', obj=ris_standard_257)


offsqs_insert_266 = [(30, seqs_30), (36, seqs_35)]
ris_insert_266 = cmso.and_pattern_match_ris(offsqs_insert_266)

h5.create_array('/', 'ris_insert_266', obj=ris_insert_266)


dg_seqs = ['TGGTGA', 'TGGCGA', 'TGGCAA', 'TGGAGA', 'TGGGGACA', 'TGATGA', 'TGGTAA', 'TGGCTA', 'TGACAA', 
          'TGACGA', 'TGCTGA', 'CGGTGA', 'CGGCAA', 'CGGCGA', 'TGCCAA',  'TGTTGA', 'TGTCAA',
          'TGTCGA', 'TGAAGA', 'TGCCGA', 'TGAGGA']
dg_offsqs = [(26, dg_seqs)]
ris_dg = cmso.and_pattern_match_ris(dg_offsqs)

h5.create_array('/', 'ris_delete_261', ris_dg)

matched_ris = np.union1d(ris_standard_257, ris_insert_266)
matched_ris = np.union1d(matched_ris, ris_dg)
ris_aligned_231_no_match_257 = np.setdiff1d(aro_w.aligned_ris, matched_ris)

h5.create_array('/', 'ris_aligned_231_no_match_257', obj=ris_aligned_231_no_match_257)


h5.close()

'''

alu offset: 48 
segment: first 
segment offset: 48
test sequence count 490,722
aligned repeat count 463,227
fraction aligned 0.94
alu offset: 216 
segment: second 
segment offset: 81
test sequence count 482,183
aligned repeat count 406,126
fraction aligned 0.84
alu offset: 139 
segment: second 
segment offset: 4
test sequence count 482,709
aligned repeat count 436,975
fraction aligned 0.91
alu offset: 231 
segment: second 
segment offset: 96
test sequence count 481,422
aligned repeat count 441,200
fraction aligned 0.92


'''











