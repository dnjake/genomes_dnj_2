# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from bokeh.models import LabelSet, ColumnDataSource
from ..gene_data import gene_data_rdr as gdr

class gene_finder_cls(object) :
    min_gene_screen_space = 60
    gene_space_dtype = np.dtype([('space_start', 'u4'), ('space_end', 'u4'), ('gene_mid_pos', 'u4'),
                                 ('gene_length', 'u4'), ('gene_symbol', 'S20')])
    def __init__(self, chrom, pos_range, screen_width, gene_y_bottom, gene_y_top) :
        self.chrom = chrom
        self.first_pos, self.last_pos = pos_range
        self.pos_length = self.last_pos - self.first_pos
        self.screen_width = screen_width
        self.gene_y_bottom = gene_y_bottom
        self.gene_y_top = gene_y_top
        self.gene_mid_y = (self.gene_y_bottom + self.gene_y_top)/2
        self.space_min_pos = self.pos_length*self.min_gene_screen_space/self.screen_width
        self.half_space_min_pos = self.space_min_pos/2
        gene_rdr = gdr.genes_rdr_cls(self.chrom)
        self.genes = gene_rdr.genes_in_interval(self.first_pos, self.last_pos)
        self.gene_spaces = np.zeros(self.genes.size, dtype=self.gene_space_dtype)
        self.gene_spaces['gene_mid_pos'] = (self.genes['tx_end'] + self.genes['tx_start'])/2
        genes_space_start = self.gene_spaces['gene_mid_pos'] - self.half_space_min_pos
        self.gene_spaces['space_start'] = genes_space_start
        self.gene_spaces['space_end'] = self.gene_spaces['gene_mid_pos'] + self.half_space_min_pos
        self.gene_spaces['gene_length'] = self.genes['tx_end'] - self.genes['tx_start']
        self.gene_spaces['gene_symbol'] = self.genes['gene_symbol']
        self.do_genes()
        
        
    def do_gene_quad_coords(self) :
        self.genes_top = None
        self.genes_bottom = None
        self.genes_left = None
        self.genes_right = None
        if self.genes is None :
            return
        self.genes_top = np.zeros(self.genes.size,dtype='i4')
        self.genes_top[:] = self.gene_y_top
        self.genes_bottom = np.zeros(self.genes.size,dtype='i4')
        self.genes_bottom[:] = self.gene_y_bottom
        self.genes_left = self.genes['tx_start']
        self.genes_right = self.genes['tx_end']        
        
        
    def assign_gene_label(self, space_start, space_end) :
        ind_start = self.gene_spaces['space_start'].searchsorted(space_start)
        ind_end = self.gene_spaces['space_end'].searchsorted(space_end)
        if ind_start == ind_end : return
        gene_spaces = self.gene_spaces[ind_start:ind_end]
        ind_choice = gene_spaces['gene_length'].argmax()
        gene_space_start, gene_space_end, gene_mid_pos, gene_length, gene_symbol = gene_spaces[ind_choice]
        symbol = str(gene_symbol)
        if len(symbol) > 10 :
            symbol = symbol[:9] + '-'
        self.label_symbols.append(symbol)
        self.label_mid_pos.append(gene_mid_pos)
        if gene_space_start - space_start > self.space_min_pos :
            self.new_spaces.append((space_start, gene_space_start))
        if space_end - gene_space_end > self.space_min_pos :
            self.new_spaces.append((gene_space_end, space_end))
        
    def do_gene_labels(self) :
        self.label_symbols = []
        self.label_mid_pos = []
        self.new_spaces = [(self.first_pos, self.last_pos)]
        while len(self.new_spaces) > 0 :
            old_spaces = self.new_spaces
            self.new_spaces = []
            for space_start, space_end in old_spaces :
                self.assign_gene_label(space_start, space_end)
        gene_count = len(self.label_symbols)
        self.label_mid_y = []
        if gene_count > 0 :
            for ind in range(gene_count) :
                self.label_mid_y.append(self.gene_mid_y)
            self.data_source = ColumnDataSource({'x': self.label_mid_pos, 'y': self.label_mid_y, 'vals': self.label_symbols})                
            self.gene_labels = LabelSet(x='x', y='y', text='vals', source=self.data_source, level='glyph',
                                        render_mode='canvas', text_baseline='middle', text_align='center',
                            text_font_size=('8pt'), text_font_style=('bold'), text_alpha=1.0)
        else :
            self.data_source = None
            self.gene_labels = None
            
    def do_genes(self) :
        self.do_gene_quad_coords()
        self.do_gene_labels()
            
            
        