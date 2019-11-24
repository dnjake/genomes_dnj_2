# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import numpy as np
import gzip
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')

'''
ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/kgXref.txt.ga
'''

'''
knownGene

name	uc001aaa.3	varchar(255)	values	Name of gene
chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
strand	+	char(1)	values	+ or - for strand
txStart	11873	int(10) unsigned	range	Transcription start position (or end position for minus strand item)
txEnd	14409	int(10) unsigned	range	Transcription end position (or start position for minus strand item)
cdsStart	11873	int(10) unsigned	range	Coding region start (or end position if for minus strand item)
cdsEnd	11873	int(10) unsigned	range	Coding region end (or start position if for minus strand item)
exonCount
exonStarts	11873,12612,13220,	longblob	 	Exon start positions (or end positions for minus strand item)
exonEnds	12227,12721,14409,	longblob	 	Exon end positions (or start positions for minus strand item)
proteinID	 	varchar(40)	values	UniProt display ID, UniProt accession, or RefSeq protein ID
alignID	uc001aaa.3	varchar(255)	values	Unique identifier (GENCODE transcript ID for GENCODE Basic)
'''

'''
kgXref

field	example	SQL type	info	description
kgID	uc001aaa.3	varchar(255)	values	Known Gene ID
mRNA	NR_046018	varchar(255)	values	mRNA ID
spID	 	varchar(255)	values	UniProt protein Accession number
spDisplayID	 	varchar(255)	values	UniProt display ID
geneSymbol	DDX11L1	varchar(255)	values	Gene Symbol
refseq	NR_046018	varchar(255)	values	RefSeq ID
protAcc	 	varchar(255)	values	NCBI protein Accession number
description	Homo sapiens DEAD/H (Asp-Gl...	longblob	 	Description
rfamAcc	 	varchar(255)	values	Rfam accession number
tRnaName	 	varchar(255)	values	Name from the tRNA track

'''


class transcript_rdr_cls(object) :
    in_transcript_dtype = np.dtype([('index', np.uint32), ('chrom', np.uint16), ('start', np.uint32), ('end', np.uint32), 
                                    ('gene_name', 'O'), ('data', 'O'), ('exons', 'O')])
    out_exon_dtype = np.dtype([('gene_region_index', np.uint32), ('transcript_index', '<u4'), ('chrom', np.uint16),
                               ('exon_start', '<u4'), ('exon_end', '<u4'), ('strand', 'S1')])    
    data_folder = 'grch37_hg19_annotation_data'
    in_file_name = 'knownGene.txt.gz'
    in_file_path = os.path.join(data_folder, in_file_name)


    def read_data(self) :    
        f = gzip.open(self.in_file_path)
        file_lines = f.readlines()
        f.close()
        transcripts = []
        for l in file_lines :
            g = l.split()
            gene_name = g[0]
            gene_symbol = ''
            chrom = g[1]
            chrom = chrom[3:]
            try:
                chrom = int(chrom)
            except ValueError:
                continue
            #index = 0
            #gene_index = 0
            strand = g[2]
            tx_start = g[3]
            tx_end = g[4]
            cds_start = g[5]
            cds_end = g[6]
            exon_count = g[7]
            exon_count = int(exon_count)
            transcript_data = (0, gene_name, gene_symbol, strand, chrom, tx_start, tx_end, cds_start, cds_end, exon_count)
            exon_starts = g[8].split(',')
            exon_starts = exon_starts[:exon_count]
            exon_ends = g[9].split(',')
            exon_ends = exon_ends[:exon_count]
            exons = []
            for s, e in zip(exon_starts, exon_ends) :
                exons.append((chrom, s, e, strand))
            exons = tuple(exons)
            transcripts.append((0, chrom, tx_start, tx_end, gene_name, transcript_data, exons ))
        self.transcripts = np.array(transcripts, dtype=self.in_transcript_dtype)

    def sort_transcripts(self) :
        self.transcripts.sort(order=['chrom', 'start', 'end'])
        indexes = np.arange(1, self.transcripts.size+1)
        self.transcripts['index'] = indexes
        names = self.transcripts['gene_name']
        gene_names = names.tolist()
        gene_name_lens = map(len, gene_names)
        self.max_gene_name_len = max(gene_name_lens)
        
    def parse_transcripts(self) :
        indexes = self.transcripts['index']
        data = self.transcripts['data']
        exons = self.transcripts['exons']
        self.full_transcripts = []
        self.full_exons = []
        for i in xrange(indexes.size) :
            index = indexes[i]
            item = data[i]
            item_exons = exons[i]
            fd = [index]
            fd.extend(item)
            self.full_transcripts.append(tuple(fd))
            for chrom, start, end, strand in item_exons :
                self.full_exons.append((0, index, chrom, start, end, strand))
                
    def generate_arrays(self, gene_name_length, gene_symbol_length) :
        gene_name_type = 'S' + str(gene_name_length)
        gene_symbol_type = 'S' + str(gene_symbol_length)
        out_transcript_dtype = np.dtype([('transcript_index', np.uint32), ('gene_region_index', np.uint32),
                            ('gene_name', gene_name_type), ('gene_symbol', gene_symbol_type), ('strand', 'S1'),
                            ('chrom', np.uint16), ('tx_start', np.uint32), ('tx_end', np.uint32),
                            ('cds_start', np.uint32), ('cds_end', np.uint32), ('exon_count', np.uint32)])
        self.transcript_zero = (0, 0, '', '', '', 0, 0, 0, 0, 0, 0)
        self.full_transcripts = np.array(self.full_transcripts, dtype=out_transcript_dtype)
        self.full_exons = np.array(self.full_exons, dtype=self.out_exon_dtype)
        
    def do_work(self) :
        self.read_data()
        self.sort_transcripts()
        self.parse_transcripts()
        
class xref_rdr_cls(object) :
    data_folder = 'grch37_hg19_annotation_data'       
    in_file_name = 'kgXref.txt.gz'
    file_path = os.path.join(data_folder, in_file_name)
    
    def read_file(self) :
        xf = gzip.open(self.file_path)
        xf_lines = xf.readlines()
        xf.close()
        name_symbols = []
        for xv in xf_lines :
            xva = xv.split('\t')
            name_symbols.append((xva[0], xva[4]))
        gene_names, gene_symbols = zip(*name_symbols)
        gene_name_lens = map(len, gene_names)
        self.max_gene_name_len = max(gene_name_lens)
        gene_name_type = 'S' + str(self.max_gene_name_len)
        gene_symbol_lens = map(len, gene_symbols)
        self.max_gene_symbol_len = max(gene_symbol_lens)
        gene_symbol_type = 'S' + str(self.max_gene_symbol_len)
        name_symbol_dtype = np.dtype([('name', gene_name_type), ('symbol', gene_symbol_type)])
        name_symbols = np.array(name_symbols, dtype=name_symbol_dtype)
        self.name_symbols = name_symbols
        
        
        
class merge_transcripts_xref_cls(object) :
    out_gene_dtype = np.dtype([('gene_region_index', np.uint32), ('chrom', np.uint16),  ('start_pos', np.uint32),
                               ('end_pos', np.uint32)])
    out_gene_zero = (0, 0, 0, 0)
    data_folder = 'grch37_hg19_annotation_data'        
    out_file_name = 'all_known_gene_transcripts.h5'
    out_file_path = os.path.join(data_folder, out_file_name)
    out_transcript_table_name = 'all_known_gene_transcripts'
    out_exon_table_name = 'transcript_exons'
    out_gene_table_name = 'known_gene_regions'
    
    def __init__(self) :
        self.trans_obj = transcript_rdr_cls()
        self.trans_obj.do_work()
        self.xref_obj = xref_rdr_cls()
        self.xref_obj.read_file()
        
        
    def set_dtypes(self) :
        gene_name_length = self.xref_obj.max_gene_name_len
        gene_symbol_length = self.xref_obj.max_gene_symbol_len
        self.trans_obj.generate_arrays(gene_name_length, gene_symbol_length)

    def set_gene_symbols(self) :
        name_symbols = self.xref_obj.name_symbols
        transcript_array = self.trans_obj.full_transcripts
        name_symbols.sort(order='name')
        transcript_array.sort(order='gene_name')
        inds_names = name_symbols['name'].searchsorted(transcript_array['gene_name'])
        m = name_symbols['name'][inds_names] == transcript_array['gene_name']
        assert m.all()
        transcript_array['gene_symbol'] = name_symbols['symbol'][inds_names]
        transcript_array.sort(order=['chrom', 'tx_start', 'tx_end'])
        self.transcript_array = transcript_array
        self.exon_array = self.trans_obj.full_exons

    def find_gene_regions(self) :
        gene_regions = [self.out_gene_zero]
        next_gene_index = 1
        for chrom in range(1,23) :
            bound = chrom + 1
            cis = self.transcript_array['chrom'].searchsorted([chrom, bound])
            chrom_transcripts = self.transcript_array[cis[0]:cis[1]]
            ind_start = 0
            ctr = chrom_transcripts[0] 
            gr_start = ctr['tx_start']
            gr_end = ctr['tx_end']
            count = 1
            for i in xrange(1, chrom_transcripts.size) :
                ctr = chrom_transcripts[i]
                tx_start = ctr['tx_start']
                tx_end = ctr['tx_end']
                if tx_start < gr_end :
                    count += 1
                    if tx_end > gr_end :
                        gr_end = tx_end
                else :
                    ind_bound = ind_start + count
                    gene_transcripts = chrom_transcripts[ind_start:ind_bound]
                    gene_region_index = next_gene_index
                    next_gene_index += 1
                    gene_regions.append((gene_region_index, chrom, gr_start, gr_end))
                    gene_transcripts['gene_region_index'] = gene_region_index
                    ind_start = i
                    gr_start = tx_start
                    gr_end = tx_end
                    count = 1
            ind_bound = ind_start + count
            gene_transcripts = chrom_transcripts[ind_start:ind_bound]
            gene_region_index = next_gene_index
            next_gene_index += 1
            gene_regions.append((gene_region_index, chrom, gr_start, gr_end))
            gene_transcripts['gene_region_index'] = gene_region_index                    
        self.gene_regions = np.array(gene_regions, dtype=self.out_gene_dtype)

    def set_exon_index_fields(self) :
        self.transcript_array.sort(order='transcript_index')
        transcript_indexes = self.transcript_array['transcript_index']
        gene_region_indexes = self.transcript_array['gene_region_index']
        self.exon_array.sort(order='transcript_index')
        tis, starts, counts = np.unique(self.exon_array['transcript_index'], return_index=True, return_counts=True)
        for ti, start, count in zip(tis, starts, counts) :
            bound = start + count
            ti_exons = self.exon_array[start:bound]
            ind_ti = transcript_indexes.searchsorted(ti)
            ti_exons['gene_region_index'] = gene_region_indexes[ind_ti]

    def write_data(self) :
        h5 = tb.open_file(self.out_file_path, 'w', filters=filters)
        out_transcript_table = h5.create_table('/', self.out_transcript_table_name,
                                               description=self.transcript_array.dtype)
        out_transcript_table.append([self.trans_obj.transcript_zero])
        out_transcript_table.append(self.transcript_array)
        out_exon_dtype = self.exon_array.dtype
        out_exon_table = h5.create_table('/', self.out_exon_table_name, description=out_exon_dtype)
        out_exon_table.append(self.exon_array)
        out_gene_table = h5.create_table('/', self.out_gene_table_name, description=self.out_gene_dtype)
        out_gene_table.append(self.gene_regions)
        h5.close()

    def do_work(self) :
        self.set_dtypes()
        self.set_gene_symbols()
        self.find_gene_regions()
        self.set_exon_index_fields()
        self.write_data()

'''        
mtxo = merge_transcripts_xref_cls()
mtxo.do_work()
'''    
    
        
        
