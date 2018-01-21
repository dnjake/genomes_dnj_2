# -*- coding: utf-8 -*-

import numpy as np
import tables as tb

filters = tb.Filters(complevel=5, complib='zlib')

gene_file_name = 'ref_seq_genes.h5'
gene_table_name = 'ref_seq_genes'

class gene_filter_cls(object) :
    def __init__(self) :
        h5 = tb.open_file(gene_file_name, 'r')
        table = getattr(h5.root, gene_table_name)
        self.gene_data = table[:]
        h5.close()

    def remove_dup_symbols(self) :        
        d = self.gene_data
        d.sort(order='gene_symbol')
        sym, index, count = np.unique(d['gene_symbol'], return_index=True, return_counts=True)
        lengths = d['tx_end'] - d['tx_start']
        data = zip(sym, index, count)
        out_data = []
        for sym, index, count in data  :
            if count == 1 :
                out_data.append(d[index])
            else :
                sym_lengths = lengths[index:index+count]
                arg_max = sym_lengths.argmax()
                max_index = index + arg_max
                out_data.append(d[max_index])                
        self.single_sym_data = np.array(out_data, d.dtype)
        self.single_sym_data.sort(order=['chrom', 'tx_start'])
    

    def filter_chrom_gene_overlaps(self, chrom_data) :
        out_data = []
        last = chrom_data[0]
        last_end = chrom_data['tx_end'][0]
        last_length = chrom_data['tx_end'][0] - chrom_data['tx_start'][0]
        for gene in chrom_data[1:] :
            chrom, tx_start, tx_end, sym = gene 
            length = tx_end - tx_start
            no_ov = False
            if tx_start > last_end :
                out_data.append(last)
                no_ov = True
            if no_ov or (length > last_length) :
                last = gene
                last_end = tx_end
                last_length = length
        out_data.append(last)
        return np.array(out_data, dtype=chrom_data.dtype)

    def filter_and_write_chrom(self, chrom) :
        bound = chrom + 1
        inds = self.single_sym_data['chrom'].searchsorted([chrom, bound])
        chrom_data =self.single_sym_data[inds[0]:inds[1]]
        filtered_data = self.filter_chrom_gene_overlaps(chrom_data)
        table_name = 'chrom_'+ str(chrom) +'_genes'
        chrom_table = self.h5.create_table('/', table_name, description=self.single_sym_data.dtype )
        chrom_table.append(filtered_data)
        chrom_table.close()
            
    def do_filter_genes(self) :
        self.remove_dup_symbols()
        self.h5 = tb.open_file('chrom_filtered_genes.h5', 'w', filters=filters)
        for chrom in range(1,23) :
            self.filter_and_write_chrom(chrom)
        self.h5.close()
            
            
            
gfo = gene_filter_cls() 
gfo.do_filter_genes()




