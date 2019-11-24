# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import gzip
import os
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')

'''
data source
ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz
'''


'''
field	example	SQL type	info	description
bin	585	smallint(6)	range	Indexing field to speed chromosome range queries.
chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
chromStart	10468	int(10) unsigned	range	Start position in chromosome
chromEnd	11240	int(10) unsigned	range	End position in chromosome
name	CpG: 115	varchar(255)	values	CpG Island
length	772	int(10) unsigned	range	Island Length
cpgNum	115	int(10) unsigned	range	Number of CpGs in island
gcNum	573	int(10) unsigned	range	Number of C and G in island
perCpg	29.8	float	range	Percentage of island that is CpG
perGc	74.2	float	range	Percentage of island that is C or G
obsExp	1.09	float	range	Ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island

'''

class cpg_writer_cls(object) :
    data_dtype = np.dtype([('index', np.uint32), ('chrom', np.uint16), ('start', np.uint32), ('end', np.uint32),
                       ('cpg_num', np.uint32), ('gc_num', np.uint32), ('per_cpg', np.float32), ('per_gc', np.float32),
                       ('obs_exp', np.float32)])
    data_folder = 'grch37_hg19_annotation_data'

    def __init__(self) :
        self.in_file_path = os.path.join(self.data_folder, self.in_file_name)

    def process_data(self) :
        f = gzip.open(self.in_file_path)
        d = f.readlines()
        f.close()
        
        out_data = [(0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0)]
        for l in d :
            lis = l.split('\t')
            chrs = lis[1]
            try :
                chrom = int(chrs[3:])
            except :
                continue
            start = int(lis[2])
            end = int(lis[3])
            cpg_num = int(lis[6])
            gc_num = int(lis[7])
            per_cpg = float(lis[8])
            per_gc = float(lis[9])
            obs_exp = float(lis[10])
            out_data.append((0, chrom, start, end, cpg_num, gc_num, per_cpg, per_gc, obs_exp))    
        out_data = np.array(out_data, dtype=self.data_dtype)
        out_data.sort(order=['chrom', 'start'])
        data_indexes = np.arange(out_data.size)
        out_data['index'] = data_indexes    
        self.out_data = out_data
        
    def write_data(self) :
        out_file_path = os.path.join(self.data_folder, self.out_file_name)
        h5 = tb.open_file(out_file_path, 'w', filters=filters)
        out_table = h5.create_table('/', self.table_name, description=self.data_dtype)
        out_table.append(self.out_data)
        h5.close()
        
    def do_work(self) :
        self.process_data()
        self.write_data()
        
class cpg_masked_writer_cls(cpg_writer_cls) :
    in_file_name = 'cpgIslandExt.txt.gz'
    out_file_name = 'cpg_masked_islands.h5'
    table_name = 'cpg_masked_islands'
    
    def __init__(self) :
        cpg_writer_cls.__init__(self)
    
    
class cpg_unmasked_writer_cls(cpg_writer_cls) :    
    in_file_name = 'cpgIslandExtUnmasked.txt.gz'
    out_file_name = 'cpg_all_islands.h5'
    table_name = 'cpg_all_islands'
    
    def __init__(self) :
        cpg_writer_cls.__init__(self)

'''
mwo = cpg_masked_writer_cls()
mwo.do_work()
awo = cpg_unmasked_writer_cls()
awo.do_work()
'''






