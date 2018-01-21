# -*- coding: utf-8 -*-

import numpy as np
from . import series_sample_data as ssd

class color_finder_cls(object) :
    def __init__(self, allele_masks) :
        self.allele_masks = allele_masks
        self.find_colors()

    def find_colors(self) :
        self.series_color = np.zeros(self.allele_masks.shape[0], dtype='S7')
        for ind in range(self.allele_masks.shape[0]) :
            self.series_color[ind] = ssd.get_color_for_allele_mask(self.allele_masks[ind])
        
            