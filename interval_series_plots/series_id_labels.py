# -*- coding: utf-8 -*-

from bokeh.models import LabelSet, ColumnDataSource

class id_labels_cls(object) :
    def __init__(self, series_snp_count, series_allele_count, series_space_start, series_y_top) :
        self.series_space_start = series_space_start
        self.label_y_bottom = series_y_top.copy()
        self.label_y_bottom += 5
        self.label_values = []
        for ind in range(series_allele_count.size) :
            self.label_values.append(str(series_snp_count[ind]) + '_' + str(series_allele_count[ind]))
        label_data = {'x': self.series_space_start, 'y': self.label_y_bottom, 'vals': self.label_values}
        self.label_source = ColumnDataSource(label_data)
        self.id_labels = LabelSet(x='x', y='y', text='vals', source=self.label_source, level='glyph', render_mode='canvas',
                          text_baseline='bottom', text_align='left',
                          text_font_size=('9pt'), text_font_style=('bold'), text_alpha=1.0)


class match_labels_cls(object) :
    def __init__(self, match_allele_count, series_space_end, series_y_top) :
        self.series_space_end = series_space_end
        self.label_y_bottom = series_y_top.copy()
        self.label_y_bottom += 5
        self.label_values = []
        for ind in range(match_allele_count.size) :
            self.label_values.append(str(match_allele_count[ind]))
        label_data = {'x': self.series_space_end, 'y': self.label_y_bottom, 'vals': self.label_values}
        self.label_source = ColumnDataSource(label_data)
        self.match_labels = LabelSet(x='x', y='y', text='vals', source=self.label_source, level='glyph', 
                                     render_mode='canvas',text_baseline='bottom', text_align='right',
                                     text_font_size=('9pt'), text_font_style=('bold'), text_alpha=1.0)
