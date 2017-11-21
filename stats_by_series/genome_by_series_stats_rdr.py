# -*- coding: utf-8 -*-

import tables as tb
import numpy as np
from . import genome_by_series_stats as bss
from . import series_stats_helpers as ssh
from ..html_display import array_table as html

series_sample_count_bins = np.array([32, 64, 128, 256, 512, 1024, 2048], dtype='u2')
series_length_bins = np.array([10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000], dtype='f4')
series_snp_count_bins = np.array([8, 16, 32, 64, 128, 256, 512], dtype='u2')



class stats_rdr_cls(object) :

    def __init__(self) :
        h5 = tb.open_file(bss.stats_file_path, 'r')
        stats_table = getattr(h5.root, bss.stats_table_name)
        self.data = stats_table[:]
        h5.close()


'''
May want to break this class into separate classes,
one for the all series table and one for the tables
of stats by pop type

Another possibility is to go back to a single table
that has pop_type, series_count, sample_count length snp count
'''

class all_series_stats_cls(object) :
    def __init__(self) :
        stats_rdr = stats_rdr_cls()
        self.data = stats_rdr.data
        self.build_stats()

    def build_stats(self) :
        stats = []
        data = self.data
        fields = ('sample_count', 'length', 'snp_count')
        for name in fields :
            field_stats = []
            field_data = data[name]
            field_stats.append(name)
            field_stats.append(np.mean(field_data))
            field_stats.append(np.median(field_data))
            field_stats.append(field_data.max())
            stats.append(tuple(field_stats))
        self.all_series_stats = stats

    def stats_html(self) :
        headers = ('stat', 'mean', 'median', 'max')
        html_table = ['<table>']
        html_table.append(html.build_simple_table_header(headers))
        for data_row in self.all_series_stats :
            table_row = ['<tr>']
            table_row.append('<td>')
            table_row.append(data_row[0])
            table_row.append('</td>')
            for d in data_row[1:] :
                table_row.append(html.num_tag)
                table_row.append(html.big_fmt.format(int(d)))
                table_row.append('</td>')
            table_row.append('</tr>')
            html_table.append(''.join(table_row))
        html_table.append('</table>')
        table_html = '\n'.join(html_table)
        return table_html       
        




class pop_type_series_stats_cls(object) :
    def __init__(self) :
        stats_rdr = stats_rdr_cls()
        self.data = stats_rdr.data
        self.build_stats()
    
    def stats_from_data(self, data, stats) :
        stats.append(data.size)
        samples = data['sample_count']
        stats.append(np.mean(samples))
        stats.append(np.median(samples))
        lengths = data['length']
        stats.append(np.mean(lengths))
        stats.append(np.median(lengths))
        stats.append(lengths.max())
        snps = data['snp_count']
        stats.append(np.mean(snps))
        stats.append(np.median(snps))
        stats.append(snps.max())

    def build_stats(self) :
        self.pop_stats = []
        stats = ['all']
        self.stats_from_data(self.data, stats)
        self.pop_stats.append(tuple(stats))
        pt = bss.pop_types
        pop_types = self.data['pop_type']
        for pop in ('no_afr', 'low_afr', 'some_afr', 'high_afr', 'all_afr') :
            t = pt[pop]
            m = pop_types == t
            pop_type_data = self.data[m]
            stats = [pop]
            self.stats_from_data(pop_type_data, stats)
            self.pop_stats.append(tuple(stats))

    def stats_html(self) :
        html_table = ['<table>']
        header_html = ['<tr>']
        header_html.append('<th rowspan="2">pop type</th>')
        header_html.append('<th rowspan="2">series count</th>')
        header_html.append('<th colspan="2" style="text-align: center;">sample count</th>')
        header_html.append('<th colspan="3" style="text-align: center;">length</th>')
        header_html.append('<th colspan="3" style="text-align: center;">snp count</th>')
        header_html.append('</tr>')
        header_html.append('<tr>')
        header_html.append('<th>mean</th>')
        header_html.append('<th>median</th>')
        header_html.append('<th>mean</th>')
        header_html.append('<th>median</th>')
        header_html.append('<th>max</th>')
        header_html.append('<th>mean</th>')
        header_html.append('<th>median</th>')
        header_html.append('<th>max</th>')
        header_html.append('</tr>')
        html_table.append(''.join(header_html))
        for dr in self.pop_stats :
            r = ['<tr>']
            r.append('<td>')
            r.append(dr[0])
            r.append('</td>')
            r.append(html.num_tag)
            r.append(html.big_fmt.format(dr[1]))
            r.append('</td>')
            for d in dr[2:4] :
                r.append(html.num_tag)
                r.append(html.int_fmt.format(int(d)))
                r.append('</td>')
            for d in dr[4:7] :
                r.append(html.num_tag)
                r.append(html.big_fmt.format(int(d)))
                r.append('</td>')
            for d in dr[7:] :
                r.append(html.num_tag)
                r.append(html.int_fmt.format(int(d)))
                r.append('</td>')
            r.append('</tr>')
            html_table.append(''.join(r))
        html_table.append('</table>')
        table_html = '\n'.join(html_table)
        return table_html       
        


class bin_stats_cls(object) :
    def __init__(self) :
        stats_rdr = stats_rdr_cls()
        self.data = stats_rdr.data
    
    def bin_counts(self, data, stats_dict) :
        length_data = data['length']
        length_bins = ssh.data_bin_counts(length_data, series_length_bins)
        stats_dict['length'].append(tuple(length_bins))
        sample_count_data = data['sample_count']
        sample_count_bins = ssh.data_bin_counts(sample_count_data, series_sample_count_bins)
        stats_dict['sample_count'].append(tuple(sample_count_bins))
        snp_count_data = data['snp_count']
        snp_count_bins = ssh.data_bin_counts(snp_count_data, series_snp_count_bins)
        stats_dict['snp_count'].append(tuple(snp_count_bins))

    def build_bin_stats(self) :
        self.bin_stats = {'pop_type':[], 'stat_max':{}, 'length':[], 'sample_count':[], 'snp_count':[]}
        length_data = self.data['length']
        self.bin_stats['stat_max']['length'] = length_data.max()
        sample_count_data = self.data['sample_count']
        self.bin_stats['stat_max']['sample_count'] = sample_count_data.max()
        snp_count_data = self.data['snp_count']
        self.bin_stats['stat_max']['snp_count'] = snp_count_data.max()
        self.bin_stats['pop_type'].append('all')
        self.bin_counts(self.data, self.bin_stats)
        pt = bss.pop_types
        pop_types = self.data['pop_type']
        for pop in ('no_afr', 'low_afr', 'some_afr', 'high_afr', 'all_afr') :
            t = pt[pop]
            m = pop_types == t
            pop_type_data = self.data[m]
            self.bin_stats['pop_type'].append(pop)
            self.bin_counts(pop_type_data, self.bin_stats)
    
    def html_table(self, pop_types, stats_data, stats_bins, max_stats_val, tag_str=None, fmt_str=None) :
        header_html = []
        header_html = ['<thead><tr><th></th>']
        stats_bins = list(stats_bins)
        stats_bins.append(max_stats_val)
        for bin_val in stats_bins :
            header_html.append('<th style="text-align:center">')
            if fmt_str is None :
                header_html.append(str(bin_val))
            else :
                header_html.append(fmt_str.format(bin_val))
            header_html.append('</th>')
        header_html.append('</tr></thead>\n')
        row_in = zip(pop_types, stats_data)
        rows_html = []
        for pop, data in row_in :
            row = ['<tr><td>', pop, '</td>']
            if tag_str is None :
                tag_str = '<td>'                        
            for di in data :
                row.append(tag_str)
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
                
    def lengths_html(self) :
        length_bins = series_length_bins.astype('u4')
        s = self.bin_stats
        length_max = s['stat_max']['length']
        num_tag = '<td style="text-align: right;">'
        big_fmt =  '{:,}'
        table_html = self.html_table(s['pop_type'], s['length'], length_bins, length_max, num_tag, big_fmt)
        return table_html
        
    def sample_count_html(self) :
        sample_count_bins = series_sample_count_bins
        s = self.bin_stats
        sample_count_max = s['stat_max']['sample_count']
        num_tag = '<td style="text-align: right;">'
        table_html = self.html_table(s['pop_type'], s['sample_count'], sample_count_bins, sample_count_max, num_tag)
        return table_html

    def snp_count_html(self) :
        snp_count_bins = series_snp_count_bins
        s = self.bin_stats
        snp_count_max = s['stat_max']['snp_count']
        num_tag = '<td style="text-align: right;">'
        table_html = self.html_table(s['pop_type'], s['snp_count'], snp_count_bins, snp_count_max, num_tag)
        return table_html



    
