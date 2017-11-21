# -*- coding: utf-8 -*-

import numpy as np
from sample_data_plus import sample_count_cls
sc = sample_count_cls()
import genome_data.sample_data_plus as sdp
import tables as tb

'''
It should be relatively easy to build allele masks for each of the plus regions
and for each country.  In fact the simplest approach might be to select allele
masks by country and then just sum them to get the one for the region.

'''

region_dtype = np.dtype([('index', 'u2'), ('code', 'S3'), ('name', 'S35'), ('allele_mask', 'u1', 5008)])

country_dtype = np.dtype([('index', 'u2'), ('region_index', 'u2'), ('code', 'S3'), 
                          ('name', 'S35'), ('allele_mask', 'u1', 5008)])


country_count = len(sdp.pop_num_map_in)
region_count = len(sdp.super_pop_plus_nums_in)

countries = np.zeros(country_count, country_dtype)
regions = np.zeros(region_count, region_dtype)

country_codes = countries['code']
country_indexes = countries['index']
country_indexes[:] = np.arange(country_count, dtype='u2')
country_region_indexes = countries['region_index']
pop_code_to_index = {}
for code, index, region_index in sdp.pop_num_map_in :
    country_codes[index] = code
    country_region_indexes[index] = region_index
    pop_code_to_index[code] = index
    
country_names = countries['name']
for code, name in sdp.pop_names.iteritems() :
    index = pop_code_to_index[code]
    country_names[index] = name
    
super_pop_code_to_index = {}
region_indexes = regions['index']
region_indexes[:] = np.arange(region_count, dtype='u2')
region_codes = regions['code']
for code, index in sdp.super_pop_plus_nums_in :
    region_codes[index] = code
    super_pop_code_to_index[code] = index
    
region_names = regions['name']
for code, name in sdp.super_pop_plus_names.iteritems() :
    index = super_pop_code_to_index[code]
    region_names[index] = name

phased_allele_calls = sc.phased_allele_calls
allele_country_indexes = phased_allele_calls['pop_num']
country_allele_masks = countries['allele_mask']
region_allele_masks = regions['allele_mask']
country_region_indexes = countries['region_index']
for country_index in range(countries.size) :
    country_allele_indexes = np.where(allele_country_indexes == country_index)[0]
    country_allele_masks[country_index, country_allele_indexes] = 1
    region_index = country_region_indexes[country_index]
    region_allele_masks[region_index, country_allele_indexes] = 1
    
    
    
filters = tb.Filters(complevel=5, complib='zlib')
out_file = 'allele_country_and_region_plus_maps.h5'
h5 = tb.open_file(out_file, 'w', filters=filters)
country_table = h5.create_table('/', 'country_alleles', description=countries.dtype)
country_table.append(countries)
region_table = h5.create_table('/', 'region_alleles', description=regions.dtype)
region_table.append(regions)
h5.close()

print 'done'
    
    
    
    
    
    
    
    
    
    