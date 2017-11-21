# -*- coding: utf-8 -*-

from chrom_snp_series_finder import chrom_snp_series_finder_cls, snp_series_table_writer_cls
chroms = range(1,23)

for chrom in chroms :
    print 'starting', chrom
    table_writer = snp_series_table_writer_cls(chrom)    
    self = chrom_snp_series_finder_cls(chrom, table_writer)
    self.organize_chrom_processing()
    table_writer.close()
    print 'completed', chrom
print 'done'