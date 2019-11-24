# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import gzip
import numpy as np
import tables as tb
filters = tb.Filters(complevel=5, complib='zlib')


data_folder = 'grch37_hg19_annotation_data'
#in_file_name = 'ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz'
#in_file_name = 'ALL.wgs.meltv1_1.20140228.ALU.low_coverage.genotypes.vcf.gz'
in_file_name = 'ALL.autosomes.pindel.20130502.complexindex.low_coverage.genotypes.vcf.gz'
in_file_path = os.path.join(data_folder, in_file_name)

f = gzip.open(in_file_path, 'r')
fl = f.readlines()
f.close()