# -*- coding: utf-8 -*-

import numpy as np
from ..autosome_snp_data.chrom_snp_series_rdr import chrom_snp_series_factory_cls
from ..html_display import array_table as html
from .. interval_series_plots import series_sample_data_table as sdt


class series_snps_cls(object) :
    unique_count_dtype = np.dtype([('count', 'u2'), ('snps', 'u2')])
    def __init__(self, chrom, series_data, series_allele_mask,
                 series_snp_data, series_snp_allele_masks) :
        self.chrom = chrom
        self.series_data = series_data
        self.series_id = str(self.series_data['item_count']) + '_' + str(self.series_data['p90_allele_count'])
        self.series_allele_mask = series_allele_mask
        self.series_snp_data = series_snp_data
        self.series_snp_allele_masks = series_snp_allele_masks
        self.not_series_allele_mask = np.logical_not(self.series_allele_mask)
        
    @classmethod
    def snps_for_series(cls, chrom, series_data, series_allele_mask ) :
        series_snps_rdr = chrom_snp_series_factory_cls(chrom)
        data = series_snps_rdr.item_objs_from_series_data(series_data)
        return cls(chrom, series_data, series_allele_mask, data.snp_data, data.allele_masks)

    def unique_snps_per_allele(self, allele_mask=None) :
        if allele_mask is None :
            allele_mask = np.ones(5008, '?')
        snp_allele_masks = np.logical_and(self.series_snp_allele_masks, allele_mask)
        snps_per_allele = snp_allele_masks.sum(axis=0)
        snps_per_allele = snps_per_allele[allele_mask]
        snp_counts, instance_counts = np.unique(snps_per_allele, return_counts=True)
        data = zip(snp_counts, instance_counts)
        data_a = np.array(data, dtype=self.unique_count_dtype)
        return data_a

    def snps_from_aps_value(self, value, allele_mask=None) :
        snp_allele_masks = self.series_snp_allele_masks
        aps = snp_allele_masks.sum(axis=0)
        m = aps == value
        if allele_mask is not None :  
            m = np.logical_and(m, allele_mask)
        snp_allele_masks = np.logical_and(snp_allele_masks, m)
        spa = snp_allele_masks.sum(axis=1)
        return spa, m

    def not_series_snp_allele_mask_from_snp_number(self, snp_number) :
        return np.logical_and(self.not_series_allele_mask, self.series_snp_allele_masks[snp_number])

    def and_not_series_allele_mask_from_snp_numbers(self, snp_numbers) :
        out_allele_mask = self.not_series_allele_mask
        for index in snp_numbers :
            index_allele_mask = self.series_snp_allele_masks[index]
            out_allele_mask = np.logical_and(out_allele_mask, index_allele_mask)
        return out_allele_mask
    
    def or_not_series_allele_mask_from_snp_numbers(self, snp_numbers) :
        first_number = snp_numbers[0]
        out_allele_mask = np.logical_and(self.not_series_allele_mask, self.series_snp_allele_masks[first_number])
        for index in snp_numbers[1:] :
            index_allele_mask = np.logical_and(self.not_series_allele_mask, self.series_snp_allele_masks[index])
            out_allele_mask = np.logical_or(out_allele_mask, index_allele_mask)
        return out_allele_mask
    
    def not_series_yes_no_allele_mask(self, yes_snp_numbers=None, no_snp_numbers=None) :
        if yes_snp_numbers is None :
            out_allele_mask = self.not_series_allele_mask
        else :
            out_allele_mask = self.and_not_series_allele_mask_from_snp_numbers(yes_snp_numbers)
        if no_snp_numbers is not None :
            out_allele_mask = np.logical_and(out_allele_mask, self.or_not_series_allele_mask_from_snp_numbers(no_snp_numbers))
        return out_allele_mask
    
    def snp_allele_mask_from_index(self, index, allele_mask=None) :
        snp_allele_masks = self.series_snp_allele_masks
        if allele_mask is not None :  
            snp_allele_masks = np.logical_and(snp_allele_masks, allele_mask)
        return snp_allele_masks[index]

    def snp_allele_counts_and_indexes(self, allele_mask) :
        snp_allele_masks = np.logical_and(self.series_snp_allele_masks, allele_mask)
        snps_per_allele = snp_allele_masks.sum(axis=0)
        m = snps_per_allele > 0
        snps_per_allele = snps_per_allele[m]
        allele_indexes = np.where(m)[0]
        return snps_per_allele, allele_indexes

    def series_alleles_per_snp(self) :
        snp_allele_masks = np.logical_and(self.series_snp_allele_masks, self.series_allele_mask)
        alleles_per_snp = snp_allele_masks.sum(axis=1)
        return alleles_per_snp
    
    def alleles_per_snp(self, allele_mask=None) :
        snp_allele_masks = self.series_snp_allele_masks
        if allele_mask is not None :
            snp_allele_masks = np.logical_and(snp_allele_masks, allele_mask)
        alleles_per_snp = snp_allele_masks.sum(axis=1)
        return alleles_per_snp

    def get_in_series_snps_per_allele(self) :
        allele_mask = self.series_allele_mask
        counts_indexes = self.snp_allele_counts_and_indexs(allele_mask)
        self.in_series_snps_per_allele, self.in_series_allele_indexes = counts_indexes

    def get_out_series_snps_per_allele(self) :
        allele_mask = self.series_allele_mask
        allele_mask = np.logical_not(allele_mask)
        counts_indexes = self.snp_allele_counts_and_indexs(allele_mask)
        self.out_series_snps_per_allele, self.out_series_allele_indexes = counts_indexes

    def snps_per_allele(self, snp_mask) :
        snp_allele_masks = self.snp_allele_masks
        snps_per_allele = snp_allele_masks[snp_mask].sum(axis=0)
        return snps_per_allele

    def series_plot_data(self) :
        series_data = np.array([self.series_data], self.series_data.dtype)
        series_allele_masks = self.series_allele_mask
        series_allele_masks.shape = (1, series_allele_masks.size)
        plot_data = {}
        for name in ['data_index', 'first_pos', 'last_pos'] :
            plot_data[name] = series_data[name]
        plot_data['series_length'] = series_data['last_pos'] - series_data['first_pos']
        plot_data['snp_count'] = series_data['item_count']
        plot_data['allele_count'] = series_allele_masks.sum(axis=1)
        plot_data['allele_mask'] = series_allele_masks
        plot_data['match_allele_count'] = None
        plot_data['match_allele_mask'] = None
        plot_data['series_data'] = series_data
        return plot_data
        
    def series_data_html(self) :
        self.plot_data = self.series_plot_data()
        self.match_test_allele_mask = None        
        data_table = sdt.series_data_table_cls(self)
        return data_table.series_data_html()

    def snp_data_html(self) :
        s_d = self.series_snp_data
        snp_data_indexes = s_d['snp_index']
        snp_numbers = np.arange(snp_data_indexes.size,dtype='i4')
        column_info = [('snp', snp_numbers, sdt.center_tag, sdt.int_fmt),
                       ('index', snp_data_indexes, sdt.num_tag, sdt.int_fmt),
                       ('pos', s_d['pos'], sdt.num_tag, sdt.big_fmt),
                       ('id', s_d['id']),
                       ('niv', s_d['not_expressed_is_variant']),  
                       ('alleles', s_d['all_count'], sdt.num_tag, sdt.int_fmt)]
        column_info.extend(sdt.region_stats_from_allele_masks(self.series_snp_allele_masks))                       
        data_html_obj = html.html_table_cls(column_info)
        html_table = data_html_obj.assemble_table()
        return html_table

    def series_and_snps_html(self) :
        series_plot_data = self.series_plot_data()
        series_snp_count = series_plot_data['snp_count'][0]
        series_sample_count = series_plot_data['allele_count'][0]
        series_id = str(series_snp_count) + '_' + str(series_sample_count)
        out_html = []
        out_html.append('<p>')
        out_html.append('<b>' + series_id + ' series</b>')
        out_html.append('</p><p>')
        out_html.append(self.series_data_html())
        out_html.append('</p>')
        out_html.append('<p>')
        out_html.append('<b>' + series_id + ' series snps</b>')
        out_html.append('</p><p>')
        out_html.append(self.snp_data_html())
        series_html = '\n'.join(out_html)
        return series_html
    
                