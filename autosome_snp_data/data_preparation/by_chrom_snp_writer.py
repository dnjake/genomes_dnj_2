# -*- coding: utf-8 -*-

# iterate through each line of a thousand genome file
# parse data
# check for validity
import os
import gzip
import numpy as np
import tables as tb
from genome_data.sample_data_plus import sample_count_cls
sc = sample_count_cls()
filters = tb.Filters(complevel=5, complib='zlib')

class thousand_genome_to_hdf5_cls(object) :
    total_allele_count = 5008
    total_sample_count = total_allele_count / 2
    min_allele_count = 16
    max_allele_count = total_allele_count - min_allele_count
    class tgd_field_indexes(object) :
        chrom = 0
        pos = 1
        id = 2
        ref = 3
        alt = 4
        qual = 5
        filter = 6
        info = 7
        format = 8
        alleles_start = 9

    snp_data_descr = [ ('snp_index', '<u4'), ('chrom', 'u2'), ('ref', 'S1'), ('alt', 'S1'), ('id', 'S11'), 
                      ('not_expressed_is_variant', 'u1'), ('pos', '<u4'), ('all_count', 'u2') ]
    snp_data_descr.extend(sc.super_pop_plus_obs_to_freq_descr)
    snp_data_descr.extend(sc.pop_count_descr)
    snp_data_dtype = np.dtype(snp_data_descr)
    snp_bitpacked_values_dtype = np.dtype([('snp_index', 'u4'), ('bitpacked_values', 'u1', (626,))])                  
    genome_data_dir = r'H:\thousand_genome_data\SNPS_1000Genome\phase3'
    genome_file_name_begin = 'ALL.chr'
    genome_file_name_end = '.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz'
    h5_file_name_begin = 'chrom_'
    h5_file_name_end = '_snp_data.h5'
    snp_data_table_name = 'snp_data'
    snp_bitpacked_allele_values_table_name = 'bitpacked_allele_values'
    def __init__(self, chrom) :
        self.chrom = chrom
        str_chrom = str(chrom)
        file_name_end = self.genome_file_name_end
        self.genome_file_name = self.genome_file_name_begin + str_chrom + file_name_end
 
        
    def is_valid_snp(self, tgd_fields) :
        tfi = self.tgd_field_indexes
        if  ((len(tgd_fields[tfi.ref])==1) and
             (len(tgd_fields[tfi.alt])==1) and
             (tgd_fields[tfi.qual]=='100') and
             (tgd_fields[tfi.filter]=='PASS')) :
                return True
        return False
        
    def get_info_fields(self, info_string) :
        vals = info_string.split(';')
        info_vals = {}
        for v in vals :
            fval = v.split('=')
            if len(fval) == 2 :
                info_vals[fval[0]] = fval[1]
        return info_vals

    def is_allele_count_valid(self, info_fields) :    
        allele_count = int(info_fields['AC'])        
        if ((allele_count < self.min_allele_count) or
            (allele_count > self.max_allele_count)):
                return False
        return True
    
    def tgd_iter(self, file_path) :
        fin = gzip.open(file_path, 'rb')
        for l in fin :
            if (l == '') or (l[0] == '#') :
                continue
            tgd_fields = l.split()
            if not self.is_valid_snp(tgd_fields) :
                continue
            info_fields = self.get_info_fields(tgd_fields[self.tgd_field_indexes.info])
            if self.is_allele_count_valid(info_fields) :            
                yield tgd_fields, info_fields
        fin.close()            
    
    def extract_allele_values(self, tgd_fields) :
        samples = tgd_fields[self.tgd_field_indexes.alleles_start:]
        assert len(samples) == self.total_sample_count
        outVals = []
        for svs in samples :
            pvs = svs.split('|')
            outVals.extend(pvs[:2])
        outA = np.array(outVals,dtype='u1')
        return outA
    
    def not_expressed_is_variant(self, allele_count) :
        if (allele_count > self.total_sample_count ) :
            return 1
        else :
            return 0

    def snp_item_data(self, snp_index, all_count, not_expressed_is_variant, tgd_fields, info_fields) :
        chrom, pos, id, ref, alt = tgd_fields[:5]
        if chrom == 'X' :
            chrom = 23
        chrom = int(chrom)
        pos = int(pos)
        data = [snp_index, chrom, ref, alt, id, not_expressed_is_variant, pos, all_count]
        return data

    def snp_stats_data(self, allele_values) :
        allele_indexes = np.where(allele_values)[0]
        super_pop_plus_obs_to_preds, pop_counts = sc.super_pop_plus_and_pop_stats(allele_indexes)
        data = list(super_pop_plus_obs_to_preds)
        data.extend(list(pop_counts))
        return data
        
    def write_data_table(self, snp_item_data, snp_data_table) :
        snp_data_table.append([snp_item_data])
        
    def write_bitpacked_allele_table(self, snp_index, allele_values, bitpacked_allele_table) :
        data_a = np.zeros(1, self.snp_bitpacked_values_dtype)
        item = data_a[0]
        item['snp_index'] = snp_index
        item['bitpacked_values'] = np.packbits(allele_values)
        bitpacked_allele_table.append(data_a)
    
    
    def write_snp_data(self, in_file_path, out_snp_data_table, out_bitpacked_allele_table) :
        for snp_index, line_data in enumerate(self.tgd_iter(in_file_path)) :
            tgd_fields, info_fields = line_data                        
            allele_count = int(info_fields['AC'])
            not_expressed_is_variant = self.not_expressed_is_variant(allele_count)
            if not_expressed_is_variant :
                all_count = self.total_allele_count - allele_count
            else :
                all_count = allele_count
            allele_values = self.extract_allele_values(tgd_fields)
            if not_expressed_is_variant :
                allele_values = np.logical_not(allele_values)
            snp_data = self.snp_item_data(snp_index, all_count, not_expressed_is_variant, tgd_fields, info_fields) 
            snp_data.extend(self.snp_stats_data(allele_values))
            self.write_data_table(snp_data, out_snp_data_table)
            self.write_bitpacked_allele_table(snp_index, allele_values, out_bitpacked_allele_table)
            if snp_index % 10000 == 0 :
                print 'writing', snp_index
            
    def do_chrom_processing(self) :
        in_file_path = os.path.join(self.genome_data_dir, self.genome_file_name)
        out_file_name = self.h5_file_name_begin + str(self.chrom) + self.h5_file_name_end
        out_h5 = tb.open_file(out_file_name, 'w', filters=filters)
        out_snp_data_table = out_h5.create_table('/', self.snp_data_table_name, self.snp_data_dtype)
        out_bitpacked_values_table = out_h5.create_table('/', self.snp_bitpacked_allele_values_table_name, 
                                                         self.snp_bitpacked_values_dtype)
        self.write_snp_data(in_file_path, out_snp_data_table, out_bitpacked_values_table)
        out_h5.close()                                                         
                    
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    