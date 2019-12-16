# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import genomes_dnj_2.repeat_analysis.offset__cluster_matches as ocm



class seq_alignment_cls(object) :
    def __init__(self, patterns, base_offsets) :
        self.patterns = patterns
        self.base_offsets = base_offsets
        self.alu_offset = base_offsets[0]
        if self.alu_offset > 135 :
            self.segment = 'second'
            self.segment_offset = self.alu_offset - 135
        else :
            self.segment = 'first'
            self.segment_offset = self.alu_offset
            
    
alignments = {}


p1 = ocm.match_pattern_cls(0, ['CGGG', 'TGGG', 'CAGG'])
p2 = ocm.match_pattern_cls(9, ['TGG'])
p3 = ocm.match_pattern_cls(12, ['CTCA'])
patterns_4 = [p1, p2, p3]
base_offsets_4 =  [4, 3, 5]
alignments[4] = seq_alignment_cls(patterns_4, base_offsets_4)    

p1 = ocm.match_pattern_cls(0, ['CGCC', 'TGCC', 'CACC'])
p2 = ocm.match_pattern_cls(6, ['TAAT', 'TAGT', 'CAAT'])
p3 = ocm.match_pattern_cls(10, ['CCCA', 'CCTA', 'CTCA', 'TCCA'])
patterns_20 = [p1, p2, p3]
base_offsets_20 =  [20, 19, 21]
alignments[20] = seq_alignment_cls(patterns_20, base_offsets_20)    

p1 = ocm.match_pattern_cls(2, ['AGG', 'AGA'])
p2 = ocm.match_pattern_cls(5, ['CAGG', 'CGGG', 'TGGG'])
p3 = ocm.match_pattern_cls(9, ['TGGAT', 'CAGAT', 'CGGAT', 'AGGAT'])
patterns_48 = (p1, p2, p3)
base_offsets_48 =  (48, 47, 49)
alignments[48] = seq_alignment_cls(patterns_48, base_offsets_48)    

p1 = ocm.match_pattern_cls(0, ['AGGTCA', 'AGCCCA', 'AGGCCA', 'AGCTCA', 'AGGTCG', 'AAGTCA', 'AGATCA'])
p2 = ocm.match_pattern_cls(6, ['GGAGTT', 'GGAGAT', 'AGAGTT', 'GGAATT', 'GAAGTT', 'AGAGAT'])
p3 = ocm.match_pattern_cls(12, ['CGAGAC', 'CAAGAC', 'TGAGAC', 'GGAGAC', 'CCAGAC'])
patterns_68 = [p1, p2, p3]
base_offsets_68 = [68, 66, 67, 69, 65]
alignments[68] = seq_alignment_cls(patterns_68, base_offsets_68)    

p1 = ocm.match_pattern_cls(2, ['AGACCA', 'AGATCA', 'AGACTA'])
p2 = ocm.match_pattern_cls(8, ['GCCTG', 'TCCTG', 'GCCTA', 'GTCTG', 'ACCTG', 'GCTTG'])
patterns_80 = [p1, p2]
base_offsets_80 =  [80, 78, 77, 79, 81]
alignments[80] = seq_alignment_cls(patterns_80, base_offsets_80)  


p1 = ocm.match_pattern_cls(0, ['AGACCA', 'AGATCA', 'AGACTA', 'AAACCA', 'AGGCCA'])
p2 = ocm.match_pattern_cls(6, ['GCCTGG', 'TCCTGG', 'GCCTGA', 'GCCTAG', 'GTCTGG', 'ACCTGG', 'GCTTGG'])
p3 = ocm.match_pattern_cls(-14, ['AGGTCAGG', 'AGCCCAGG', 'AGGCCAGG', 'AGGTCAAG',
                               'AGCTCAGG', 'AGGTCGGG', 'AAGTCAGG', 'AGATCAGG', 'AGGTCAGA'])
patterns_82 = [p1, p2, p3]
base_offsets_82 =  [82, 80, 77, 81, 83, 79]
alignments[82] = seq_alignment_cls(patterns_82, base_offsets_82)  

'''
Out[19]: 'CAACATGGTGAAACCC'

       (1094377493,  91,   523, 78988), (1094377493,  92,  1642, 78988),
       (1094377493,  93, 17448, 78988), (1094377493,  94,  4624, 78988),
       (1094377493,  95, 46738, 78988), (1094377493,  96,  2329, 78988),
       (1094377493,  97,   668, 78988), (1094377493,  98,   428, 78988),

'''

'''
Out[22]: 'GGTGAAACCCCATCTC'

       (2919322845,  80,    14, 40796), (2919322845,  81,    10, 40796),
       (2919322845,  82,    34, 40796), (2919322845,  83,    54, 40796),
       (2919322845,  84,    32, 40796), (2919322845,  85,    48, 40796),
       (2919322845,  86,    52, 40796), (2919322845,  87,   166, 40796),
       (2919322845,  88,    47, 40796), (2919322845,  89,    64, 40796),
       (2919322845,  90,    90, 40796), (2919322845,  91,   125, 40796),
       (2919322845,  92,   141, 40796), (2919322845,  93,   162, 40796),
       (2919322845,  94,   137, 40796), (2919322845,  95,   204, 40796),
       (2919322845,  96,   267, 40796), (2919322845,  97,   305, 40796),
       (2919322845,  98,  1016, 40796), (2919322845,  99, 13857, 40796),
       (2919322845, 100,  2174, 40796), (2919322845, 101, 19364, 40796),
       (2919322845, 102,  1009, 40796), (2919322845, 103,   314, 40796),
       (2919322845, 104,   207, 40796), (2919322845, 105,   214, 40796),
       (2919322845, 106,    98, 40796), (2919322845, 107,    80, 40796),
       (2919322845, 108,    56, 40796), (2919322845, 109,    55, 40796),
       (2919322845, 110,    27, 40796), (2919322845, 111,    39, 40796),
 
'''

p1 = ocm.match_pattern_cls(0, ['GGTGAA', 'AGTGAA', 'GGCAAA', 'GGAGAA', 'GGCGAA', 'AGTGAG', 'GATGAA', 'GGTAAA', 'AGCAAG'])
p2 = ocm.match_pattern_cls(2, ['TGAAACCC', 'CAAAACCC', 'AGAAACCC', 'CGAAACCC', 'TGAAACCT', 'TGAAACTC', 'TGAGACCC', 
                               'TGAAATCC', 'TAAAACCC', 'CAAGACCC'])
p3 = ocm.match_pattern_cls(4, ['AAACCCCG', 'AAACCCCA', 'AAACCCTG', 'AGACCCTG', 'AGACCCCA', 'AAACCTCA', 'AAACCTCG', 
                                'AAACTCCA', 'AAATCCCA', 'AAACTCTG'])
p4 = ocm.match_pattern_cls(6, ['ACCCCATC', 'ACCCCGTC', 'ACCCTGTC', 'ACCTCATC', 'ACCTCGTC',  'ACTCCATC', 'ATCCCATC', 
                               'ACCCTATC',  'ACCTTGTC', 'ACTCTGTC'])
patterns_101 = [p1, p2, p3, p4]
base_offsets_101 =  [101, 99, 100, 102, 98, 97, 103, 96, 104, 105, 95, 94, 93, 92, 91, 90, 105, 106, 107, 108, 109, 110]
alignments[101] = seq_alignment_cls(patterns_101, base_offsets_101)    


p1 = ocm.match_pattern_cls(0, ['TACTAAAAATACAAAA'])
p2 = ocm.match_pattern_cls(-6, ['TGTCTC', 'CATCTC', 'CGTCTC', 'TATCTC'])
p3 = ocm.match_pattern_cls(-12, ['AAACCCCA', 'AAACCCTG', 'AAACCCCG', 'AGACCCTG', 'AGACCCCA', 
                                  'AAACCTCA', 'AAACTCCA','AAACCTCG', 'AAACCCTA',  'AAATCCCA', 'AAACCTTG'])
patterns_117 = [p1, p2, p3]
base_offsets_117 = [117, 115, 116, 114, 118, 113, 119, 112, 120, 121, 111, 110, 108, 107, 106, 105, 104, 103,
               122, 123, 124, 125]

alignments[117] = seq_alignment_cls(patterns_117, base_offsets_117)    
    

p1 = ocm.match_pattern_cls(-5, ['TTAGCCG', 'TTAGCCA', 'TTAGCTG', 'TTAGCTA', 'TTAGCCT'])
p2 = ocm.match_pattern_cls(6, ['TGGTGG', 'TAGTGG', 'CGGTGG', 'TGGTAG', 'TGATGG', 'TGGTGA',
                               'CAGTGG', 'TGTTGG', 'TTGTGG', 'TGGTGT', 'TGGTTG'])
p3 = ocm.match_pattern_cls(2, ['GGCATGG','GGCGTGG', 'GGTGTGG', 'GGTATGG', 'GATGTGG',
                                'AGCATGG', 'GGCATGA', 'AGTGTGG', 'GGCTTGG', 'GGCATAG', 'GACATGG'])
patterns_139 = [p1, p2, p3]
base_offsets_139 = [139, 138, 137, 140, 141, 142, 136, 143, 135, 134, 143, 144, 133, 131, 
                    130, 129, 128, 127, 126, 125, 124, 123, 122, 121, 145, 146, 147, 148, 149]
alignments[139] = seq_alignment_cls(patterns_139, base_offsets_139)    

p1 = ocm.match_pattern_cls(0, ['TGTAAT', 'TGTAGT', 'TGTGGT', 'TATAAT', 'TATAGT'])
p2 = ocm.match_pattern_cls(6, ['CCCAGC', 'CCTAGC', 'CCCAGT', 'CCCAAC', 'CTCAGC', 'TCCAGC'])
p3 = ocm.match_pattern_cls(10, ['GCTACT', 'GCTATT', 'GCTGCT'])
patterns_159 = [p1, p2, p3]
base_offsets_159 = [159, 158, 157, 160, 161, 162, 156, 163, 155, 163, 164, 154, 153, 152, 151, 165, 166, 167,
               151, 150, 149, 148, 147, 146, 168, 169]
alignments[159] = seq_alignment_cls(patterns_159, base_offsets_159)    


p1 = ocm.match_pattern_cls(-3, ['CGGGAGGC', 'TGGGAGGC', 'CAGGAGGC'])
p2 = ocm.match_pattern_cls(5, ['TGAGGCA', 'TGAGGTG', 'TGAGACA', 'TGAAGCA', 'TAAGGCA', 'TGAGGTA', 'TGAGGCG', 'CGAGGCA'])
p3 = ocm.match_pattern_cls(10, ['CAGGAGAAT', 'TGGGAGGAT', 'CAGGAGGAT', 'CAGAAGAAT', 'CAAGAGAAT',
                                'CATGAGAAT', 'CAGGAGAAC', 'CACAAGAAT', 'CAGGAAAAT', 'CACGAGAAT'])
patterns_178 = [p1, p2, p3]
base_offsets_178 = [178, 177, 176, 179, 180, 175, 181, 174, 182, 173, 183, 184, 185, 186, 187, 
                165, 172, 171, 170, 169, 168, 188, 167, 166, 189]
alignments[178] = seq_alignment_cls(patterns_178, base_offsets_178)    



p1 = ocm.match_pattern_cls(0, ['GAGGTTG', 'GAGCTTG', 'GAGGCTG', 'GAAGTTG', 'GAGATTG'])
p2 = ocm.match_pattern_cls(-15, ['TGAACCCA', 'TGAACCTG', 'TGAACCCG', 'TGAGCCCA', 'TGAGCCTG'])
p3 = ocm.match_pattern_cls(-9, ['CAGGAGGCA', 'CGGGAGGCG', 'TGGGAGGCA', 'CAGGAGGCG', 'TGGGAGGCG',
                                'CGGGAGGCA', 'TGGGAGGTG', 'CAGGAGGTG', 'CGGGAGGTG'])
p4 = ocm.match_pattern_cls(7, ['CAGTGAGC', 'CGGTGAGC', 'CAGTGAGT', 'CAGTAAGC', 'CAGTGAAC', 'TAGTGAGC', 'CAATGAGC'])

patterns_216 = [p1, p2, p3, p4]
base_offsets_216 = [216, 215, 214, 217, 218, 213, 219, 212, 220, 211, 221, 210, 222, 209, 223, 224]
alignments[216] = seq_alignment_cls(patterns_216, base_offsets_216)    



p1 = ocm.match_pattern_cls(0, ['TGCACT', 'TGTACT', 'TGCATT', 'TACACT', 'CGCACT'])
p2 = ocm.match_pattern_cls(6, ['CCAGCC', 'CTAGCC', 'CCAGTC', 'CCAGCT', 'CCAACC', 'TCAGCC', 'GCAGCC'])
patterns_245 = [p1, p2]
base_offsets_245 = [245, 244, 243, 246, 247, 248, 242, 241, 249, 250, 240, 239, 251, 252, 
                253, 239, 238, 237, 236, 235, 234, 233, 232, 254, 255, 256, 257, 258, 231, 230]
alignments[245] = seq_alignment_cls(patterns_245, base_offsets_245)    

p1 = ocm.match_pattern_cls(0, ['CGAGAT', 'TGAGAT', 'CAAGAT', 'TATGAT', 'CATGAT',  'TGTGAT', 'AGAGAT', 'GGAGAT'])
p2 = ocm.match_pattern_cls(-4, ['GAGCCG', 'GAGCCA', 'GAGCTG', 'GAGCTA', 'GAGCAG'])
p3 = ocm.match_pattern_cls(2, ['AGATCA', 'AGATCG', 'AGATTG', 'TGATCA', 'AGATGG', 'TGATTG', 
                                'AGATAG', 'AGATCC', 'AGATCT', 'TGATCG', 'AGATTA'])
p4 = ocm.match_pattern_cls(8, ['TGCCAC',  'CACCAC', 'CGCCAC', 'TGCCAT', 'CGCCAT', 'CACCAT', 'TACCAC'])
patterns_231 = [p1, p2, p3, p4]
base_offsets_231 = [offset-14 for offset in base_offsets_245]
alignments[231] = seq_alignment_cls(patterns_231, base_offsets_231)


'''
n = centers['center'][65]

Out[14]: 'ACTCCGTCTCAAAAAA'

       (493735936, 254,   72, 35577), (493735936, 255,   69, 35577),
       (493735936, 256,   83, 35577), (493735936, 257,  104, 35577),
       (493735936, 258,  116, 35577), (493735936, 259,  140, 35577),
       (493735936, 260,  712, 35577), (493735936, 261,  211, 35577),
       (493735936, 262,  248, 35577), (493735936, 263,  242, 35577),
       (493735936, 264,  250, 35577), (493735936, 265,  261, 35577),
       (493735936, 266,  278, 35577), (493735936, 267,  361, 35577),
       (493735936, 268,  461, 35577), (493735936, 269,  776, 35577),
       (493735936, 270, 1965, 35577), (493735936, 271, 3462, 35577),
       (493735936, 272, 8119, 35577), (493735936, 273, 6073, 35577),
       (493735936, 274, 2997, 35577), (493735936, 275, 1794, 35577),
       (493735936, 276, 1109, 35577), (493735936, 277,  737, 35577),
       (493735936, 278,  522, 35577), (493735936, 279, 1380, 35577),
       (493735936, 280,  431, 35577), (493735936, 281,  336, 35577),
       (493735936, 282,  279, 35577), (493735936, 283,  205, 35577),
       (493735936, 284,  183, 35577), (493735936, 285,  118, 35577),
       (493735936, 286,  105, 35577), (493735936, 287,   91, 35577),
       (493735936, 288,   86, 35577), (493735936, 289,   92, 35577),
       (493735936, 290,   63, 35577), (493735936, 291,   66, 35577),


'''


'''
n = centers['center'][14]
Out[13]: 'GTCTCAAAAAAAAAAA'

       (3074424832, 250,    78, 81483), (3074424832, 251,   100, 81483),
       (3074424832, 252,    86, 81483), (3074424832, 253,   126, 81483),
       (3074424832, 254,   128, 81483), (3074424832, 255,   108, 81483),
       (3074424832, 256,   128, 81483), (3074424832, 257,   139, 81483),
       (3074424832, 258,   175, 81483), (3074424832, 259,   223, 81483),
       (3074424832, 260,   236, 81483), (3074424832, 261,   266, 81483),
       (3074424832, 262,   314, 81483), (3074424832, 263,   348, 81483),
       (3074424832, 264,   420, 81483), (3074424832, 265,  1087, 81483),
       (3074424832, 266,   606, 81483), (3074424832, 267,   608, 81483),
       (3074424832, 268,   719, 81483), (3074424832, 269,   687, 81483),
       (3074424832, 270,   786, 81483), (3074424832, 271,   886, 81483),
       (3074424832, 272,  1089, 81483), (3074424832, 273,  1406, 81483),
       (3074424832, 274,  2199, 81483), (3074424832, 275,  4525, 81483),
       (3074424832, 276,  7637, 81483), (3074424832, 277, 13922, 81483),
       (3074424832, 278, 14965, 81483), (3074424832, 279,  7311, 81483),
       (3074424832, 280,  4493, 81483), (3074424832, 281,  2719, 81483),
       (3074424832, 282,  1921, 81483), (3074424832, 283,  1350, 81483),
       (3074424832, 284,  2124, 81483), (3074424832, 285,  1042, 81483),
       (3074424832, 286,   852, 81483), (3074424832, 287,   729, 81483),
       (3074424832, 288,   572, 81483), (3074424832, 289,   469, 81483),
       (3074424832, 290,   358, 81483), (3074424832, 291,   278, 81483),
       (3074424832, 292,   279, 81483), (3074424832, 293,   250, 81483),
       (3074424832, 294,   223, 81483), (3074424832, 295,   178, 81483),
       (3074424832, 296,   166, 81483), (3074424832, 297,   139, 81483),
       (3074424832, 298,   116, 81483), (3074424832, 299,    99, 81483),

'''


p1 = ocm.match_pattern_cls(-1, ['TGTCTCAA', 'CGTCTCAA', 'CATCTCAA', 'AGTCTCAA', 'TATCTCAA', 'TGTCTCCA','TGTCTAAA', 'GGTCTCAA', 
                                'TGTCTAAA', 'AGTCTCAA', 'TGCCTCAA', 'TGTCTCTA', 'TGTCTCAG', 'TGTCTTAA', 'CATCTAAA', 'GGTCTCAA',
                                'TGTCAAAA', 'TGTTTCAA', 'CATCTCAG', 'CATCTCCA', 'CATCTTAA', 'CGTCTAAA', 'CATTTCAA', 'CACCTCAA',
                                'TGTCTCTT', 'CGTCTCAG', 'CGTCTCCA', 'TGTCTCGA', 'TGTCTGAA'])
p2 = ocm.match_pattern_cls(-4, ['CTCTGTC', 'CTCCATC', 'CTCCGTC', 'CCCTGTC''GACTGTC', 'CTCTATC', 'TTCTGTC', 'CCTTGTC', 
                                'CTCAGTC', 'CTTTGTC', 'CTCTGCC', 'TTCCATC', 'CTTCATC', 'CCCTATC', 'TCCTGTC',
                                'CTCCATT', 'CTCTGTT', 'CCCCATC', 'CTCCACC', 'CTTCGTC', 'CCCTGCC', 'CTGTGTC', 
                                'TTCCGTC', 'CCCTGTT', 'CACTGTC', 'CTCCGTT', 'CTCGGTC', 'CTCCGCC', 'CTCCCTC', 
                                'ACCTGTC', 'CTCCTTC'])
p3 = ocm.match_pattern_cls(-8, ['GAGACTCTG', 'GAGACTCCA', 'GAGACTCCG', 'AAGACTCTG', 'GAGACCCTG', 'AAGACTCCA', 'AAGACTCCG', 
                                'GAAACTCCA', 'GAAACTCTG', 'GAAACTCCG',  'AAGACCCTG', 'GAGACTCTA', 'AAAACTCTG', 'GAGACTCAG',
                                'GAGACCCTA', 'AAAACTCCG', 'GAGACCCCA', 'AAGACTCTA', 'GAGACTCCC', 'GAGACTCCT', 'GAGATTCTG',
                                'GTGAGACTG', 'GAGACTTTG', 'GCAAGACTG', 'GAGACCTTG', 'AAGACTCAG', 'GAGACTCTC', 'GAGACTCTT'])
patterns_278 = [p1, p2, p3]
base_offsets_278 = [278, 277, 276, 279, 275, 280, 281, 274, 282, 283, 284, 273, 272, 265, 285, 
                    271, 286, 270, 287, 288, 289, 290, 291, 292, 293, 294, 269, 268, 267, 266]
alignments[278] = seq_alignment_cls(patterns_278, base_offsets_278)


p1 = ocm.match_pattern_cls(-3, ['TCTCAAAA', 'TCTCCAAA', 'CCTCAAAA', 'TCTAAAAA', 'TCTCTAAA', 'TCTCAGAA', 'TCTCAAAC',
                                'TCTTAAAA', 'TTTCAAAA', 'TCTCAAGA', 'TCTCAAAG', 'TCTCAAAT'])
p2 = ocm.match_pattern_cls(-5, ['TGTCTC', 'CATCTC', 'CGTCTC', 'TATCTC', 'AGTCTC','TGCCTC', 'GGTCTC'])
p3 = ocm.match_pattern_cls(-11, ['AGACTCTG', 'AGACTCCA', 'AGACTCCG', 'AGACCCTG', 'AAACTCCA', 'AAACTCTG', 'AAACTCCG',
                                 'AGACTCTA', 'AGACTTTG', 'AGACCCTA', 'AGACTCAG', 'AGACCTTG', 'AGATTCTG',  'AGACTTCA', 
                                 'AGACCCCA', 'AGATCCTG', 'AGATTCCA'])
patterns_282 = [p1, p2, p3]
base_offsets_282 = [282, 281, 280, 283, 279, 284, 278, 285, 286, 287, 277, 288, 289,
                276, 275, 290, 274, 291, 273, 272, 271, 270, 269, 292, 293]
alignments[282] = seq_alignment_cls(patterns_282, base_offsets_282)    



































    
    