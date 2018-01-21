# -*- coding: utf-8 -*-

int_fmt =  '{:d}'
big_fmt =  '{:,}'
float_fmt = '{:.2f}'
exp_fmt = '{:.1e}'

plot_width = 1600

bottom_margin_height = 50
def bottom_margin_div() :
    return '<div style="min-height:' +str(bottom_margin_height) + 'px;"></div>'

font_size = 'font-size:10pt;'

border_collapse = 'border-collapse:collapse;'

top_vertical_align = 'vertical-align:text-top;'

fixed_table_layout = 'table-layout:fixed;'

right_align = 'text-align:right'
center_align = 'text-align:center'

table_tag = '<table style="' + font_size + border_collapse + top_vertical_align + '">'

right_align_cell_tag = '<td style="' + right_align + '">'
center_align_cell_tag = '<td style="' + center_align + '">'

