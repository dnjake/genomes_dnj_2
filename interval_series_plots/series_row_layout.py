# -*- coding: utf-8 -*-


'''
Approach will be to layout series in order of allele count.  Try to use one row
as long as possible.  But, once layout has moved on to the next row, don't go
back.  Keep track of max alleles in a series for each row.  That number will be
useful when processing gets to height determination.  Alternatively could figure
out the height for each series first and keep the max height for each row.

An array with the row for each series has to be kept in series data index order.
Think I want both the snp + allele_count id and the match count on the top
Need minimum width for series to allow for both labels.  Probably should center
short series in minimum width.  Probably means I want a start pos and an end pos
for the series space that can be used when it is time to plot the labels
'''
from __future__ import division
import numpy as np

class row_layout_cls(object) :
    series_layout_dtype = np.dtype([('row_index', 'u2'), ('pos_start', 'u4'), ('pos_end', 'u4')])
    row_layout_dtype = np.dtype([('row_index', 'u2'), ('height_allele_count', 'u2'), ('series_spaces', 'O')])
    with_match_labels_min_series_screen_space = 100
    only_series_labels_min_series_screen_space = 50
    def __init__(self, series_first_pos, series_last_pos, series_height_allele_count, series_order_allele_count,
                 x_axis_pos_range, x_axis_screen_width=900, has_matches=True) :
        self.series_first_pos = series_first_pos
        self.series_last_pos = series_last_pos
        self.series_height_allele_count = series_height_allele_count
        self.series_order_allele_count = series_order_allele_count
        self.x_pos_start, self.x_pos_end = x_axis_pos_range
        self.pos_range = float(self.x_pos_end - self.x_pos_start)
        self.x_axis_screen_width = x_axis_screen_width
        self.has_matches = has_matches
        self.min_series_screen_space = self.with_match_labels_min_series_screen_space
        if not has_matches :
            self.min_series_screen_space = self.only_series_labels_min_series_screen_space
        self.series_min_pos = float(self.min_series_screen_space*self.pos_range)/float(x_axis_screen_width)
        self.series_min_pos = int(self.series_min_pos)
        self.half_series_min_pos = self.series_min_pos/2
        self.do_layout()
        
    def find_series_row_space(self, series_array_index) :
        series_length = self.series_lengths[series_array_index]
        if series_length >= self.series_min_pos :
            space_start = self.series_first_pos[series_array_index]
            space_end = self.series_last_pos[series_array_index]
        else :
            space_start = self.series_mid_pos[series_array_index] - self.half_series_min_pos
            if space_start < self.x_pos_start:
                space_start = self.x_pos_start
            space_end = space_start + self.series_min_pos
            if space_end > self.x_pos_end:
                space_end = self.x_pos_end
                space_start = space_end - self.series_min_pos
        series_row_space = 0, space_start, space_end
        return series_row_space
    
    def layout_series(self, series_array_index) :
        series_row_space = self.find_series_row_space(series_array_index)
        series_row, series_start, series_end = series_row_space
        series_allele_count = self.series_height_allele_count[series_array_index]
        row_index, row_allele_count, row_spaces = self.layout_row
        series_row_space = row_index, series_start, series_end
        for ind, space in enumerate(row_spaces) :
            space_start, space_end = space[1:]
            '''
            The assumption at this point is that series_start is
            after the last series end
            '''
            if series_end < space_start :
                '''
                Insert series before the current enumerated space
                The series row index is already set to the right value
                '''
                self.series_row_layout[series_array_index] = series_row_space
                if ind == 0 :
                    new_row_spaces = [series_row_space]
                else :
                    new_row_spaces = row_spaces[:ind]
                    new_row_spaces.append(series_row_space)
                new_row_spaces.extend(row_spaces[ind:])
                if series_allele_count > row_allele_count :
                    row_allele_count = series_allele_count
                self.layout_row = (row_index, row_allele_count, new_row_spaces)
                return
            elif series_start < space_end :
                '''
                The new series does not end before the start of the current space
                and does start befor the end of that space.  Therefore the new
                series space overlaps the current series space and a new row must
                be started.
                '''
                self.layout_rows.append(self.layout_row)
                row_index += 1
                row_allele_count = series_allele_count
                series_row_space = row_index, series_start, series_end
                self.layout_row = (row_index, row_allele_count, [series_row_space])
                self.series_row_layout[series_array_index] = series_row_space
                return
            '''
            If processing gets here the new series does not preceed or overlap the 
            current series space.  Processing moves on to the next series
            '''
        row_spaces.append(series_row_space)
        self.series_row_layout[series_array_index] = series_row_space
        if series_allele_count > row_allele_count :
            row_allele_count = series_allele_count
            self.layout_row = (row_index, row_allele_count, row_spaces)
                
    def do_layout(self) :
        self.layout_rows = []
        self.series_mid_pos = (self.series_first_pos + self.series_last_pos)/2
        self.series_lengths = self.series_last_pos - self.series_first_pos
        self.series_row_layout = np.zeros(self.series_first_pos.size, self.series_layout_dtype)
        self.row = self.series_row_layout['row_index']
        self.pos_start = self.series_row_layout['pos_start']
        self.pos_end = self.series_row_layout['pos_end']
        self.layout_series_order = self.series_order_allele_count.argsort()
        self.layout_series_order = self.layout_series_order[::-1]
        series_array_index = self.layout_series_order[0]
        series_row_space = self.find_series_row_space(series_array_index)
        series_allele_count = self.series_height_allele_count[series_array_index]
        self.layout_row = (0, series_allele_count, [series_row_space])
        self.series_row_layout[series_array_index] = series_row_space
        for index in self.layout_series_order[1:] :
            self.layout_series(index)
        self.layout_rows.append(self.layout_row)
        self.layout_row = None
        self.layout_rows = np.array(self.layout_rows, dtype=self.row_layout_dtype)
        
        
                 