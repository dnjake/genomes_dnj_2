# -*- coding: utf-8 -*-

import tables as tb
import numpy as np
from . import genome_all_snp_stats as ass
from ..html_display import array_table as html



class stats_rdr_cls(object) :
    data_dtype = np.dtype([('min_snps_in_group', 'u2'), ('group_count', 'u4'), ('snp_count', 'u4')])    
    def __init__(self) :
        h5 = tb.open_file(ass.stats_file_path, 'r')
        stats_table = getattr(h5.root, ass.stats_table_name)
        self.data = stats_table[:]
        h5.close()

    def do_snp_stats(self) :
        d = self.data
        snp_stats = []
        snp_1_series = d['snp_1'].sum()
        snp_1_count = snp_1_series
        snp_stats.append(('1 snp', snp_1_series, snp_1_count))
        snp_2_series = d['snp_2'].sum()
        snp_2_count = 2*snp_2_series
        snp_stats.append(('2 snps', snp_2_series, snp_2_count))
        snp_3_series = d['snp_3'].sum()
        snp_3_count = 3*snp_3_series
        snp_stats.append(('3 snps', snp_3_series, snp_3_count))
        snp_3_or_less_series = snp_1_series + snp_2_series + snp_3_series
        snp_3_or_less_count = snp_1_count + snp_2_count + snp_3_count
        snp_stats.append(('3 snps or less', snp_3_or_less_series, snp_3_or_less_count ))
        snp_4_or_more_series = d['series_count'].sum()
        snp_4_or_more_count = d['series_snps'].sum()
        snp_stats.append(('4 snps or more', snp_4_or_more_series, snp_4_or_more_count))
        self.snp_stats = snp_stats
        self.snp_stats_table_header = ('series class', 'series count',  'series snps count')
        total_snps = snp_3_or_less_count + snp_4_or_more_count
        self.total_snps_table_header = ('analysis total snps', 'in series of 4 or more snps', 'series count')
        self.total_snps_data = (total_snps, snp_4_or_more_count, snp_4_or_more_series)
        
    def total_stats_html(self) :
        html_table = ['<table>']
        html_table.append(html.build_simple_table_header(self.total_snps_table_header))
        table_row = ['<tr>']
        for d in self.total_snps_data :
            table_row.append(html.center_tag)
            table_row.append(html.big_fmt.format(int(d)))
            table_row.append('</td>')
        html_table.append(''.join(table_row))
        html_table.append('</table>')
        table_html = '\n'.join(html_table)
        return table_html       
        
    def by_snp_stats_html(self) :
        html_table = ['<table>']
        html_table.append(html.build_simple_table_header(self.snp_stats_table_header))
        for r in self.snp_stats :
            row = ['<tr>']
            row.extend(('<td>',r[0],'</td>'))
            for d in r[1:] :
                row.extend((html.num_tag, html.big_fmt.format(int(d)), '</td>'))
            row.append('</tr>')
            html_table.append(''.join(row))
        html_table.append('</table>')
        table_html = '\n'.join(html_table)
        return table_html

