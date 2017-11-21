# -*- coding: utf-8 -*-

from autosome_snp_data.chrom_snp_series_rdr import chrom_snp_series_tables_cls
from chrom_snp_series_data_writer import  chrom_snp_series_data_writer_cls

chroms = range(1,23)
for chrom in chroms :
    print 'starting', chrom
    series_rdr = chrom_snp_series_tables_cls(chrom)
    data_writer = chrom_snp_series_data_writer_cls(chrom, series_rdr)
    data_writer.do_write()
    series_rdr.close()
    print 'completed', chrom
print 'done'
