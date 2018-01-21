# -*- coding: utf-8 -*-
import numpy as np
import tables as tb
import gzip

'''
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

Database: hg19    Primary Table: refGene    Row Count: 69,723   Data last updated: 2017-11-21
Format description: A gene prediction with some additional info.
field    example    SQL type    info    description
bin    585    smallint(5) unsigned    range    Indexing field to speed chromosome range queries.
name    NR_148357    varchar(255)    values    Name of gene (usually transcript_id from GTF)
chrom    chr1    varchar(255)    values    Reference sequence chromosome or scaffold
strand    +    char(1)    values    + or - for strand
txStart    11868    int(10) unsigned    range    Transcription start position (or end position for minus strand item)
txEnd    14362    int(10) unsigned    range    Transcription end position (or start position for minus strand item)
cdsStart    14362    int(10) unsigned    range    Coding region start (or end position for minus strand item)
cdsEnd    14362    int(10) unsigned    range    Coding region end (or start position for minus strand item)
exonCount    3    int(10) unsigned    range    Number of exons
exonStarts    11868,12612,13220,    longblob         Exon start positions (or end positions for minus strand item)
exonEnds    12227,12721,14362,    longblob         Exon end positions (or start positions for minus strand item)
score    0    int(11)    range    score
name2    LOC102725121    varchar(255)    values    Alternate name (e.g. gene_id from GTF)
cdsStartStat    unk    enum('none', 'unk', 'incmpl', 'cmpl')    values    enum('none','unk','incmpl','cmpl')
cdsEndStat    unk    enum('none', 'unk', 'incmpl', 'cmpl')    values    enum('none','unk','incmpl','cmpl')
exonFrames    -1,-1,-1,    longblob         Exon frame {0,1,2}, or -1 if no frame for exon
'''

filters = tb.Filters(complevel=5, complib='zlib')

genes_dtype = np.dtype([('chrom', 'u2'), ('tx_start', 'u4'), ('tx_end', 'u4'), ('gene_symbol', 'S20')])

class field_indexes(object) :
    bin = 0
    name = 1
    chrom = 2
    strand = 3
    tx_start = 4
    tx_end = 5
    cds_start = 6
    cds_end = 7
    exon_count = 8
    exon_starts = 9
    exon_ends = 10
    score = 11
    symbol = 12
    cds_start_status = 13
    cds_end_status = 14
    exon_frames = 15
fi = field_indexes
file_path = 'refGene.txt.gz'
f = gzip.open(file_path, 'rb')
genes = []
for l in f :
    fields = l.split()
    assert len(fields) == 16
    chrom = fields[fi.chrom]
    try :
        chrom_int = int(chrom[3:])
    except :
        pass
    else :
        gene_sym = fields[fi.symbol][:20]
        data = (chrom_int, int(fields[fi.tx_start]), int(fields[fi.tx_end]), gene_sym)
        genes.append(data)
f.close()

genes = np.array(genes,dtype=genes_dtype)
genes.sort(order=['chrom', 'tx_start'])

h5 = tb.open_file('ref_seq_genes.h5', 'w', filters=filters)
table = h5.create_table('/', 'ref_seq_genes', description=genes.dtype)
table.append(genes)
h5.close()