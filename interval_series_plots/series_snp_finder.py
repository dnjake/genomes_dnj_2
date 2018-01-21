# -*- coding: utf-8 -*-

from ..autosome_snp_data.chrom_snp_series_rdr import chrom_snp_series_factory_cls


class snp_finder_cls(object) :
    def __init__(self, chrom, series_data) :
        self.chrom = chrom
        self.series_data = series_data
        self.find_snps()
        
    def find_snps(self) :
        series_snps_rdr = chrom_snp_series_factory_cls(self.chrom)
        series_snps_pos = []
        for item_data in self.series_data :
            item_snps = series_snps_rdr.item_objs_from_series_data(item_data)
            series_snps_pos.append(item_snps.snp_data['pos'])
        series_snps_rdr.close()
        self.series_snps_pos = series_snps_pos

        
class snp_line_plot_cls(object) :
    def __init__(self, series_snps_pos, series_y_bottom, series_y_top) :
        self.snp_x_vals = []
        self.snp_y_vals = []
        self.colors = 'black'
        for ind, item_snp_pos in enumerate(series_snps_pos) :
            item_y_val = (series_y_bottom[ind], series_y_top[ind])
            #item_snp_pos = series_snps_pos[ind]
            for item_pos in item_snp_pos :
                item_pos_x_val = (item_pos, item_pos)
                self.snp_x_vals.append(item_pos_x_val)
                self.snp_y_vals.append(item_y_val)
        
        