# -*- coding: utf-8 -*-

import genomes_dnj_2.stats_by_pos.genome_by_series_first_last_stats as fls

gw = fls.genome_series_event_stats_writer_cls()
gw.do_all_chrom_stats()
