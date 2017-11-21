# -*- coding: utf-8 -*-

import numpy as np
from IPython.display import HTML
from bokeh.plotting import output_notebook, show
from bokeh.layouts import column, layout
from bokeh.models import Div
from ..chrom_plots import chrom_plots as cp
from ..genome_plots import genome_plot as gp
from ..stats_by_series import genome_by_series_stats_rdr as genome_series_stats
from ..series_plots import series_plot
from ..stats_by_series import genome_length_snp_series_stats as length_snp_stats

from ..stats_by_series import genome_all_snp_stats_rdr as snp_stats_rdr
snp_stats = snp_stats_rdr.stats_rdr_cls()
snp_stats.do_snp_stats()
html = ['<div style="width:700px;">']
html.append(snp_stats.total_stats_html())
html.append(snp_stats.by_snp_stats_html())
html.append('</div>')
snp_stats_html = '\n'.join(html)

from ..stats_by_series import genome_by_series_stats_rdr as series_stats_rdr
    
series_stats = series_stats_rdr.all_series_stats_cls()
html = ['<div style="width:700px;">']
html.append(series_stats.stats_html())
html.append('</div>')
all_series_stats_html = '\n'.join(html)

from ..genome_plots import genome_table_stats as gts
gsh = gts.genome_stats_html_cls()
html = ['<div style="width:700px;">']
html.append(gsh.genome_stats_html())
html.append('</div>')
genome_stats_html = '\n'.join(html)

snp_stats_plot = layout([ 
                        [Div(text=snp_stats_html)]
                      ], 
                      sizing_mode='fixed')

all_series_stats_plot = layout([ 
                        [Div(text=all_series_stats_html)]
                        ], 
                        sizing_mode='fixed')

genome_stats_html_plot = layout([ 
                            [Div(text=genome_stats_html)]
                            ], 
                            sizing_mode='fixed')

def genome_stats_plots() :
    plt0 = gp.do_plot('samples_in_series')
    plt1 = gp.do_log_plot('series_count')
    plt2 = gp.do_log_plot('sample_weighted_length')
    plt3 = gp.do_log_plot('sample_weighted_snps')
    plt4 = gp.do_log_plot('snps_in_series')
    plt5 = gp.do_log_plot('mean_series_length')
    plt6 = gp.do_log_plot('mean_series_snps')
    return column(plt0, plt1,plt2,plt3,plt4,plt5,plt6)


chrom_1 = 1
cp1 = cp.chrom_plot_cls(chrom_1)

chrom_11 = 11
cp11 = cp.chrom_plot_cls(chrom_11)


pop_type_series_stats = series_stats_rdr.pop_type_series_stats_cls()
html = ['<div style="width:700px;">']
html.append(pop_type_series_stats.stats_html())
html.append('</div>')
pop_type_series_stats_html = '\n'.join(html)

pop_type_series_stats_plot = layout([ 
                                        [Div(text=pop_type_series_stats_html)]
                                    ], 
                                    sizing_mode='fixed')

def html_h4(htext) :
    return '<h4>' + htext + '</h4>'

def pop_type_bin_stats_html() :
    bs = genome_series_stats.bin_stats_cls()
    bs.build_bin_stats()    
    html = ['<div style="width:700px;">']
    html.append(html_h4('length'))
    html.append(bs.lengths_html())
    html.append('<p>')
    html.append(html_h4('sample_count'))
    html.append(bs.sample_count_html())
    html.append('<p>')
    html.append(html_h4('snp_count'))
    html.append(bs.snp_count_html())
    html.append('</div>')
    return '\n'.join(html)


pop_type_bin_stats_plot = layout([ 
                                    [Div(text=pop_type_bin_stats_html())]
                                 ], 
                                 sizing_mode='fixed')

    
spo = series_plot.data_plot_cls()

def genome_length_snp_series_stats_html() :
    lss = length_snp_stats.length_snp_bin_stats_cls()
    html = ['<div style="width:700px;">']
    html.append(lss.html_stats())
    html.append('</div>')
    return '\n'.join(html)

genome_length_snp_series_stats_plot = layout([ 
                                                  [Div(text=genome_length_snp_series_stats_html())]
                                              ], 
                                              sizing_mode='fixed')
    
    