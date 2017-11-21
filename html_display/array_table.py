# -*- coding: utf-8 -*-
'''
The simplest approach might be just to build the data column by column
zip the columns
and iterater over the rows
'''
#import numpy as np
from itertools import izip

num_tag = '<td style="text-align: right;">'
center_tag = '<td style="text-align: center;">'
int_fmt =  '{:d}'
big_fmt =  '{:,}'
float_fmt = '{:.2f}'
exp_fmt = '{:.1e}'

def int_cells(row_html, data) :
    for d in data :
        row_html.append(num_tag)
        row_html.append(int_fmt.format(d))
        row_html.append('</td>')
        
def big_int_cells(row_html, data) :
    for d in data:
        row_html.append(num_tag)
        row_html.append(big_fmt.format(d))
        row_html.append('</td>')

def float_cells(row_html, data) :
    for d in data:
        row_html.append(num_tag)
        row_html.append(float_fmt.format(d))
        row_html.append('</td>')


def build_html_table_column(column_data, html_tag=None, format_str=None) :            
    out_html = []
    if html_tag is None :
        html_tag = '<td>'
    for item in column_data :
        if format_str is None :
            cell_content = str(item)
        else :
            cell_content = format_str.format(item)
        cell_str = html_tag + cell_content + '</td>'
        out_html.append(cell_str)
    return tuple(out_html)


def build_html_table_rows(columns_html) :
    html_rows = []
    for ind, row_cells in enumerate(izip(*columns_html)) :
        if ind % 2 == 0 :
            row_tag = '<tr class="even">'
        else :
            row_tag = '<tr>'
        content_html = ''.join(row_cells) 
        row_html = row_tag + content_html + '</tr>'
        html_rows.append(row_html)
    return html_rows

'''
could just put in a test for a tuple here
'''
def build_simple_table_header(column_names) :
    header_fields = ['<thead>', '<tr>']
    for name in column_names :
        tag = '<th style="text-align:center">'
        if type(name) is tuple :
            name, spanv = name
            tag = '<th style="text-align:center" colspan=' + str(spanv) + '>'
        header_fields.append( tag + name + '</th>')
    header_fields.extend(['</tr>', '</thead>'])
    header_html = ''.join(header_fields)
    return header_html

    
class html_table_cls(object) :
    def __init__(self, column_info) :
        self.column_info = column_info
        self.header_html = None
        self.table_rows_html = None
        
    def build_table(self) :
        column_html = []
        column_names = []
        for column in self.column_info :
            column_names.append(column[0])
            # could deal with multiple cols here
            if type(column[1]) is tuple :
                for col_data in column[1] :
                    column_html.append(build_html_table_column(*col_data))
            else :
                column_html.append(build_html_table_column(*column[1:]))
        self.header_html = build_simple_table_header(column_names)
        self.table_rows_html = build_html_table_rows(column_html)
            
    def headers_table(self) :
        if self.header_html is None :
            self.build_table()
        html_table = ['<table>']
        html_table.append(self.header_html)
        html_table.append('</table>')
        self.table_html = '\n'.join(html_table)
        return self.table_html                    

    def assemble_table(self) :
        if self.table_rows_html is None :
            self.build_table()
        html_table = ['<table>']
        html_table.append(self.header_html)
        html_table.extend(self.table_rows_html)
        html_table.append('</table>')
        self.table_html = '\n'.join(html_table)
        return self.table_html                    

    
'''
Basically what I want is a class for building an html table from an ndarray or maybe
even a dictionary that is driven by a description of the fields to include and formatting
information for each data item
'''

'''
class html_table(object) :
    class column_format(object) :
        def __init__(self, format_str, html_str=None) :
            self.format_str = format_str
            self.html_str = html_str
    integer_format = '{:d}'
    float_format = '{:.2f}'
    string_format = '{:s}'
    pad_left_format = column_format(string_format, '<td class=pad_left>')
    
    data_col_descr = namedtuple('data_col_descr',['field_name', 'field_header', 'field_format'])
    multi_col_descr = namedtuple('multi_col_descr',['multi_name', 'multi_header', 'content_descrs', 
                                          'content_data_columns', 'content_header_rows'])
    data_col = namedtuple('data_col', ['col_number', 'data_array', 'col_html', 'col_format'])
    header_cell = namedtuple('header_cell',['header_text', 'row_number', 'col_span', 'row_span' ])
    snp_stats_data_col_descrs = [
                     data_col_descr('first_snp_index', 'index', integer_format), 
                     data_col_descr('first_snp_pos', 'start', integer_format),
                     data_col_descr('series_pos_length', 'length', integer_format),
                     data_col_descr('series_snp_count', 'snps', integer_format),
                     data_col_descr('variant_density', 'density', float_format),
                     data_col_descr('all_sample_count', 'all haplotypes', integer_format),
                     data_col_descr('all_spec_mean', 'all mean', float_format),
                     data_col_descr('all_spec_min', 'all min', float_format),
                     data_col_descr('p90_sample_count', '90% haplotypes', integer_format),
                     data_col_descr('p90_spec_mean', '90% mean', float_format),
                     data_col_descr('p90_spec_min', '90% min', float_format)
                 ]
    pop_data_cols = [
            data_col_descr('haplotype_count', 'haplo', integer_format), 
            data_col_descr('haplotype_freq', 'freq', float_format), 
            data_col_descr('variant_homozygotes', 'homo', integer_format), 
            data_col_descr('heterozygotes', 'hetero', integer_format)
    ]                     
    single_pop_data_cols = [
            data_col_descr('first_snp_index', 'index', integer_format), 
            data_col_descr('series_snp_count', 'variants', integer_format)
    ]      
    pop_codes = ['ALL', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    pop_titles = ['All', 'Africa', 'America', 'East Asia', 'Europe', 'South Asia']
    pop_fields = zip(pop_codes, pop_titles)
    pop_descrs = [multi_col_descr(key, header, pop_data_cols, 4, 2) for key, header in pop_fields]
    pop_stats_descrs = []
    pop_stats_descrs.extend(single_pop_data_cols)
    pop_stats_descrs.extend(pop_descrs)
    snp_details_descrs = [
        data_col_descr('snp_index', 'index', integer_format),
        data_col_descr('pos', 'position', integer_format),
        data_col_descr('id', 'id', pad_left_format) 
    ]
    
    #table_style = ''
    #    td {padding-left: 10;
    #        padding-right: 10; }    
    #    tr.even {background-color: #efefef}
    #    ''
    def __init__(self, data, table_description) :
        self.data = data
        self.table_description = table_description
        self.table_style = table_description['html_stylesheet']
        self.data_cols = None
        self.header_cells = None
        self.data_col_count = 0
        self.data_header_row_count = 0
        self.total_header_rows = None
        self.html_header_lines = None
        self.html_data_lines = None
        self.html_table_lines = None
    
    def initialize_structures(self) :        
        self.data_cols = []
        self.header_cells = {}
        self.find_header_cell_counts()
        self.total_header_rows = self.data_header_row_count
        table_header_title = None
        first_data_row = 0
        if 'table_header_title' in self.table_description :
            self.total_header_rows += 1
            header_text = self.table_description['table_header_title']
            row_number = 0
            col_span = self.data_col_count
            row_span = 1
            table_header_title = self.header_cell(header_text, row_number, col_span, row_span)
        for ind in range(self.total_header_rows) : 
            self.header_cells[ind] = []
        if table_header_title is not None :
            self.header_cells[0].append(table_header_title)
            first_data_row = 1
        descrs = self.table_description['data_descrs']
        self.find_data_cols(descrs, self.data)
        self.find_header_cells(descrs, first_data_row)

    def find_header_cell_counts(self) :
        descrs = self.table_description['data_descrs']
        for descr in descrs :
            if type(descr) is self.data_col_descr :
                self.data_col_count += 1
                if self.data_header_row_count < 1 :
                    self.data_header_row_count = 1
            else :
                self.data_col_count += descr.content_data_columns
                if descr.content_header_rows > self.data_header_row_count :
                    self.data_header_row_count = descr.content_header_rows
    # I want to add the html string to this tupple                
    def find_data_cols(self, descrs, data) :
        for descr in descrs :
            if type(descr) is self.data_col_descr :
                col_number = len(self.data_cols)
                col_data = data[descr.field_name]
                field_format = descr.field_format
                if type(field_format) is str :
                    col_html = '<td align="right">'
                    col_format = field_format
                else :
                    col_html = field_format.html_str
                    #col_html = '<td class=' + col_css_class + '>'
                    col_format = field_format.format_str
                self.data_cols.append(self.data_col(col_number, col_data, col_html, col_format))
            else :
                field_name = descr.multi_name
                multi_data = data[field_name]
                self.find_data_cols(descr.content_descrs, multi_data)
                
    def find_header_cells(self, descrs, row_number) :
        for descr in descrs :
            if type(descr) is self.data_col_descr :
                header_text = descr.field_header
                col_span = 1
                row_span = self.total_header_rows - row_number
                self.header_cells[row_number].append(self.header_cell(header_text, row_number, col_span, row_span))
            else :
                header_text = descr.multi_header
                col_span = descr.content_data_columns
                row_span = 1
                self.header_cells[row_number].append(self.header_cell(header_text, row_number, col_span, row_span))
                self.find_header_cells(descr.content_descrs, row_number+1)

    def build_html_headers(self) :
        self.html_header_lines = []
        self.html_header_lines.append('<thead>')
        for ind in range(self.total_header_rows) :
            header_line = ['<tr>']
            for cell in self.header_cells[ind] :
                cell_str = '<th'
                if cell.col_span > 1 :
                    cell_str = cell_str + ' colspan=' + str(cell.col_span)
                if cell.row_span > 1 :
                    cell_str = cell_str + ' rowspan=' + str(cell.row_span)
                cell_str = cell_str + '>' + cell.header_text + '</th>'
                header_line.append(cell_str)
            header_line.append('</tr>')
            self.html_header_lines.append(''.join(header_line))
        self.html_header_lines.append('</thead>')

    def build_data_lines(self, ind_start, ind_bound) :
        self.html_data_lines = []
        for ind in xrange(ind_start, ind_bound) :
            data_line = []
            if ind % 2 == 0 :
                data_line.append('<tr class="even">')
            else :
                data_line.append('<tr>')
            # I probably want to make the cell string something that can be specified in the format    
            for dc in self.data_cols :
                #cell_str = '<td align="right">'
                #cell_str = '<td>'
                cell_str = dc.col_html
                cell_str = cell_str + dc.col_format.format(dc.data_array[ind]) + '</td>'
                data_line.append(cell_str)
            data_line.append('</tr>')
            self.html_data_lines.append(''.join(data_line))

    def build_html_table(self) :
        if self.data_cols is None :
            self.initialize_structures()
        if self.html_header_lines is None :
            self.build_html_headers()
        if self.html_data_lines is None :
            self.build_data_lines(0, self.data.size)
        self.html_table_lines = ['<table cellspacing=3 width="100%">']
        self.html_table_lines.extend(self.html_header_lines)
        self.html_table_lines.extend(self.html_data_lines)
        self.html_table_lines.append('</table>')
        if 'table_legend_html' in self.table_description :
            self.html_table_lines.append(self.table_description['table_legend_html'])
        
    def table_html(self) :
        if self.html_table_lines is None :
            self.build_html_table()
        return ''.join(self.html_table_lines)

    def print_table(self, table_file) :
        if self.html_table_lines is None :
            self.build_html_table()
        doc = gui.QTextDocument()
        doc.setDefaultStyleSheet(self.table_style)
        doc.setHtml(self.table_html())
        pixel_width = self.table_description['pixel_width']
        printer = gui.QPrinter()
        printer.setOutputFileName(table_file)
        printer.setPaperSize(core.QSizeF(pixel_width, 800.0), gui.QPrinter.DevicePixel)
        printer.setPageMargins(20.0, 20.0, 20.0, 30.0, gui.QPrinter.DevicePixel)
        doc.print_(printer)

    '''



