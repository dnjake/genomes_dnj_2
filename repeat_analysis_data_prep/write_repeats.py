# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import gzip
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')

'''
data source
ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
'''



'''
field	example	SQL type	description
bin	585	smallint(5) unsigned	Indexing field to speed chromosome range queries.
swScore	1504	int(10) unsigned	Smith Waterman alignment score
milliDiv	13	int(10) unsigned	Base mismatches in parts per thousand
milliDel	4	int(10) unsigned	Bases deleted in parts per thousand
milliIns	13	int(10) unsigned	Bases inserted in parts per thousand
genoName	chr1	varchar(255)	Genomic sequence name
genoStart	10000	int(10) unsigned	Start in genomic sequence
genoEnd	10468	int(10) unsigned	End in genomic sequence
genoLeft	-249240153	int(11)	-#bases after match in genomic sequence
strand	+	char(1)	Relative orientation + or -
repName	(CCCTAA)n	varchar(255)	Name of repeat
repClass	Simple_repeat	varchar(255)	Class of repeat
repFamily	Simple_repeat	varchar(255)	Family of repeat
repStart	1	int(11)	Start (if strand is +) or -#bases after match (if strand is -) in repeat sequence
repEnd	463	int(11)	End in repeat sequence
repLeft	0	int(11)	-#bases after match (if strand is +) or start (if strand is -) in repeat sequence
id	1	char(1)	First digit of id field in RepeatMasker .out file. Best ignored.


'''



class field_offsets(object) :
    sw_score = 1
    mismatches = 2
    deleted = 3
    inserted = 4
    chrom = 5
    start = 6
    end = 7
    strand = 9
    repeat_name = 10
    repeat_class = 11
    repeat_family = 12
    repeat_start = 13
    repeat_end = 14
    repeat_left = 15


class rmsk_data_writer_cls(object) :
    data_dtype = np.dtype([('index', np.uint32), ('chrom', np.uint16), ('start_pos', np.uint32),
               ('end_pos', np.uint32), ('strand', 'S1'), ('repeat_name', 'S20' ),
               ('repeat_class', 'S20'), ('repeat_family', 'S20'), ('sw_score', np.uint16),
               ('mismatches', np.uint16), ('deleted', np.uint16), ('inserted', np.uint16),
               ('repeat_start', np.int16), ('repeat_end', np.int16), ('repeat_left', np.int16)])
    data_folder = 'grch37_hg19_annotation_data'
    in_file_name = 'rmsk.txt.gz'
    in_file_path = os.path.join(data_folder, in_file_name)
    write_file_name = 'indexed_sorted_repeats.h5'
    out_file_path = os.path.join(data_folder, write_file_name)
    write_table_name = 'indexed_sorted_repeats'
    null_repeat = (0,0,0,0,'','','','',0,0,0,0,0,0,0)
    
    def read_rmsk_data(self) :
        out_data = [self.null_repeat]
        index = 0
        f = gzip.open(self.in_file_path, 'r')        
        for l in f :
            dil = l.split('\t')
            chrom = dil[field_offsets.chrom]
            try :
                chrom = int(chrom[3:])
            except :
                continue
            d_start = int(dil[field_offsets.start])
            d_end = int(dil[field_offsets.end])
            start_pos = d_start
            if start_pos < d_end :
                end_pos = d_end
            else :
                start_pos = d_end
                end_pos = d_start
            strand = dil[field_offsets.strand]
            repeat_name = dil[field_offsets.repeat_name]
            repeat_name = repeat_name[:20].lower()
            repeat_class = dil[field_offsets.repeat_class]
            repeat_class = repeat_class[:20].lower()
            repeat_family = dil[field_offsets.repeat_family]
            repeat_family = repeat_family[:20].lower()
            sw_score = int(dil[field_offsets.sw_score])
            mismatches = int(dil[field_offsets.mismatches])
            deleted = int(dil[field_offsets.deleted])
            inserted = int(dil[field_offsets.inserted])
            repeat_start = int(dil[field_offsets.repeat_start])
            repeat_end = int(dil[field_offsets.repeat_end])
            repeat_left = int(dil[field_offsets.repeat_left])
            oi = (index, chrom, start_pos, end_pos, strand, repeat_name, repeat_class, repeat_family,
                  sw_score, mismatches, deleted, inserted, repeat_start, repeat_end, repeat_left)
            out_data.append(oi)            
        f.close()
        self.repeats = np.array(out_data, dtype=self.data_dtype)
        self.repeats.sort(order=['chrom', 'start_pos'] )
        repeat_index = np.arange(self.repeats.size)
        self.repeats['index'] = repeat_index
        
            
    def write_output(self) :
        h5 = tb.open_file(self.out_file_path, 'w', filters=filters)
        out_table = h5.create_table('/', self.write_table_name, description=self.repeats.dtype)
        out_table.append(self.repeats)
        h5.close()
        
    def do_work(self) :
        self.read_rmsk_data()
        self.write_output()

'''        
rwo = rmsk_data_writer_cls()
rwo.do_work()
'''
