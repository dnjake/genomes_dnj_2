# -*- coding: utf-8 -*-

import numpy as np
#from . import genome_by_series_stats as bss
from . import series_stats_helpers as ssh
from . import genome_by_series_stats_rdr as ssr

'''
This class is for the table that correlates
length and snp distributions.
'''
class length_snp_bin_stats_cls(object) :
    length_bins = ssr.series_length_bins
    snp_count_bins = ssr.series_snp_count_bins
    def __init__(self) :
        stats_rdr = ssr.stats_rdr_cls()
        self.data = stats_rdr.data
        self.max_length = self.data['length'].max()
        self.max_snp_count = self.data['snp_count'].max()
        
    def find_length_bins(self) :
        self.data.sort(order='length')
        bin_starts = self.data['length'].searchsorted(self.length_bins)
        bin_indexes = [0]
        bin_indexes.extend(bin_starts)
        bin_indexes.append(self.data.size)
        bin_indexes = np.array(bin_indexes, dtype='i4')
        self.bin_data = np.zeros(bin_indexes.size-1, dtype='O')
        for ind in range(bin_indexes.size-1) :
            start = bin_indexes[ind]
            bound = bin_indexes[ind+1]
            self.bin_data[ind] = self.data[start:bound]
            
    def find_snp_bin_stats(self) :
        self.bin_stats = np.zeros(self.bin_data.size, dtype='O')
        for ind in range(self.bin_data.size) :
            snp_count_data = self.bin_data[ind]['snp_count']
            self.bin_stats[ind] = ssh.data_bin_counts(snp_count_data, self.snp_count_bins)

    def html_table(self) :
        length_bins = list(self.length_bins.astype('u4'))
        length_bins.append(self.max_length)
        num_tag = '<td style="text-align: right;">'
        snp_count_bins = list(self.snp_count_bins)
        snp_count_bins.append(self.max_snp_count)
        header_html = []
        header_html = ['<thead><tr><th style="text-align:center"></th>']
        for bin_val in snp_count_bins :
            header_html.append('<th style="text-align:center">')
            header_html.append(str(bin_val))
            header_html.append('</th>')
        header_html.append('</tr></thead>\n')
        row_in = zip(length_bins, self.bin_stats)
        rows_html = []
        for length_val, data in row_in :
            row = ['<tr>']
            row.append(num_tag)
            row.append('{:,}'.format(length_val))
            row.append('</td>')
            for di in data :
                row.append(num_tag)
                row.append('{:,}'.format(di))
                row.append('</td>')
            row.append('</tr>\n')
            rows_html.append(''.join(row))
        html_table = ['<table>']
        html_table.extend(header_html)
        html_table.extend(rows_html)
        html_table.append('</table>')
        table_html = ''.join(html_table)
        return table_html

    def html_stats(self) :
        self.find_length_bins()
        self.find_snp_bin_stats()
        return self.html_table()