# -*- coding: utf-8 -*-

import numpy as np

class y_values_finder_cls(object) :
    series_y_value_dtype = np.dtype([('y_bottom', 'f4'), ('y_top', 'f4')])
    top_margin = 25
    bottom_margin = 25
    row_margin = 10
    def __init__(self, row_height, series_row, series_height, label_height=25) :
        self.row_height = row_height
        self.series_row = series_row
        self.series_height = series_height
        self.label_height = label_height
        self.find_y_values()
        
    def find_row_y_values(self) :
        row_y_bottom = self.row_height.copy()
        height_delta = self.label_height + self.row_margin
        row_y_bottom += height_delta
        y_value = self.bottom_margin
        rev_row_y_bot = row_y_bottom[::-1]
        for ind in range(rev_row_y_bot.size) :
            h = rev_row_y_bot[ind]
            rev_row_y_bot[ind] = y_value
            y_value += h
        self.row_y_bottom = row_y_bottom
        self.plot_top = y_value + self.top_margin
        
    def find_series_y_values(self) :
        self.series_y_values = np.zeros(self.series_height.size, dtype=self.series_y_value_dtype)
        for ind in range(self.series_row.size) :
            series_row_index = self.series_row[ind]
            series_y_bottom = self.row_y_bottom[series_row_index]
            series_y_top = self.series_height[ind] + series_y_bottom
            self.series_y_values[ind] = (series_y_bottom, series_y_top)
            
    def find_y_values(self) :
        self.find_row_y_values()
        self.find_series_y_values()
        