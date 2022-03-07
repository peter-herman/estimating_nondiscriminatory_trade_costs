__Author__ = "Peter Herman"
__Project__ = "A Pragmatic Approach to Estimating Nondiscriminatory Non-tariff Trade Costs"
__Created__ = "August 19, 2021"
__Description__ = '''This script prepares the data for analysis.'''


import pandas as pd
import re
import zipfile as zipfile
import os as os
from economic_analysis_tools.data_analysis.data_diagnostics import CompareIdentifiers
import numpy as np


# ---
# Specification
# ----

# Gravity variables to include
grav_vars = ['distance', 'contiguity', 'common_language', 'agree_pta', 'agree_fta', 'agree_cu',
             'colony_of_destination_ever', 'colony_of_origin_ever', 'member_eu_joint', 'member_wto_joint']

# Save location
save_local = "analysis/constructed_data"

# Years to include
years = list(range(2010,2016))

# Paths for MFN tariff data
zipped_mfn_local = "data\source_data\TRAINS Bulk MFN 2005-2019"
extracted_mfn_local = os.path.join(zipped_mfn_local,"Extracted")
mfn_country_codes_local = "data\source_data\TRAINS Bulk MFN 2005-2019\\00_TRAINS_country_codes.csv"

# Path for the FAO to HS concordance
fao_hs_cncd_local = "data\source_data\\FAOSTAT_expanded_item_codes_8-19-2021.csv"

# Path for info on EUN membership
eun_members_local = "data\source_data\EUN_member_codes.csv"

# DGD gravity files
grav_local = ["data/source_data/dgd_2.1/release_2.1_2005_2009.csv",
              "data/source_data/dgd_2.1/release_2.1_2010_2014.csv",
              "data/source_data/dgd_2.1/release_2.1_2015_2019.csv"]

# NTM Data
ntm_local ='data\source_data/ntm_hs6_2010_2018v.12/NTM_hs6_2010_2018 v.12.csv'

# GDP data
gdp_local = "data\source_data\WDI GDP\WDI GDP data 08-05-2021.csv"

# ITPD bilateral trade data
itpd_local = "data\source_data\itpd_e_r01\ITPD_E_R01.csv"
itpd_concord_local = "data\source_data\itpd_e_r01\ITPD Agriculture FCL concordance.txt"



# --------------------------
# Prep Unilateral Data
# --------------------------

# ----
# HS concordance to ITPD sectors
# ----

# Load ITPD codes
itpd_fcl_raw = pd.read_csv(itpd_concord_local)



#++
# Step 1: FAO Item code to HS concordance
#++

# Get the list of FAO Item Codes in Trade data and convert to string
use_code_list = itpd_fcl_raw['fcl'].unique().tolist()
use_code_list = [str(code) for code in use_code_list]

# Load concordance
fcl_hs_ccrd_raw = pd.read_csv(fao_hs_cncd_local, dtype=str)
fcl_hs_cncd = fcl_hs_ccrd_raw.loc[fcl_hs_ccrd_raw['Item Code'].isin(use_code_list), :].copy()

# Split multiple HS codes into different rows
hs_code_dfs = list()
for col in ['HS07 Code', 'HS12 Code']:
    # Grab Item code and correspodning HS revision codes
    fao_hs_nomen = fcl_hs_cncd[['Item Code', col]].copy()
    # Expand HS codes into new columns
    hs_expanded = fao_hs_nomen[col].str.split(', ', expand = True)
    # Add expanded HS codes back on to Item codes
    hs_wide = pd.concat([fao_hs_nomen['Item Code'], hs_expanded], axis = 1)
    # Reshape HS codes into long format
    hs_wide.set_index('Item Code', inplace = True)
    hs_long = hs_wide.stack()
    hs_long = hs_long.reset_index()
    # Drop unwanted info from restack and rename columns
    hs_long.drop(['level_1'], axis = 1, inplace = True)
    hs_long.rename(columns = {0:'hs6'}, inplace = True)
    # Create column specifying nomenclature revision
    hs_long['nomen'] = col[0:4]
    # Store data in list
    hs_code_dfs.append(hs_long)

# Combine all nomenclatures
fao_hs = pd.concat(hs_code_dfs, axis = 0)
# Drop duplicates, which occurs due to the source FAO concordance listing identical FCL/HS codes under multiple
#   different "Domains"
fao_hs.drop_duplicates(inplace=True)
fao_hs.sort_values(['Item Code','hs6'], inplace = True)

# Count number of HS6 codes in each item code by HS revision
fao_hs['hs6_count'] = fao_hs.groupby(['Item Code','nomen'])['hs6'].transform('count')

# List of all needed HS6 codes
all_hs = fao_hs['hs6'].unique().tolist()

#++
# Map FCL codes to ITPD sectors
#++
itpd_fcl = itpd_fcl_raw.copy()
itpd_fcl['fcl']= itpd_fcl['fcl'].astype(str)
itpd_hs_full = itpd_fcl.merge(fao_hs, how ='left', left_on = ['fcl'], right_on = ['Item Code'],
                              validate = '1:m', indicator=True)
# Drop the one Item Code without any corresponding HS codes (code 654, "Dregs from brewing, distillation")
itpd_hs_full = itpd_hs_full.loc[itpd_hs_full['_merge'] == 'both', :]
itpd_hs_full.drop(['_merge'], axis = 1, inplace = True)

# One HS code may map to more than one FCL code, then back to one ITPD code, causing ITPD--HS duplicates.
#   Need to drop duplicates from the FCL step
itpd_hs = itpd_hs_full[['itpd', 'itpd_description', 'hs6', 'nomen']].copy()
itpd_hs.drop_duplicates(inplace = True)

itpd_hs['hs6_count'] = itpd_hs.groupby(['itpd','nomen'])['hs6'].transform('count')
itpd_hs['itpd'] = itpd_hs['itpd'].astype(str)

# ---
# MFN Tariff Data
# ---



##
# Load MFN tariff data
##

#Unzip the zipped MFN data files (note that the 2019 files exist but contain no data other than headers)
zipped_mfn = os.listdir(zipped_mfn_local)
zipped_mfn = [zp for zp in zipped_mfn if zp.endswith('.zip')]

for folder in zipped_mfn:
    with zipfile.ZipFile(os.path.join(zipped_mfn_local,folder),"r") as zip_ref:
        zip_ref.extractall(extracted_mfn_local)


# Load MFN data files
mfn_files = os.listdir(extracted_mfn_local)
mfn_files = [file for file in mfn_files if file.endswith(".CSV")]
all_sheets = list()
for num, file in enumerate(mfn_files):
    if num%300 == 0:
        print("Loading MFN: {}% complete".format(round(100*num/len(mfn_files))))
    single_sheet = pd.read_csv(os.path.join(extracted_mfn_local, file), dtype={'ProductCode':str})
    all_sheets.append(single_sheet)
mfn_raw = pd.concat(all_sheets, axis = 0)

# Drop un-needed hs codes
mfn_sub = mfn_raw.loc[mfn_raw['ProductCode'].isin(all_hs),:]


##
# Map HS MFN info to FAO Item Codes
##

# Define HS nomenclature code matching (FAO only has H2/HS07 and H3/HS12 so others are mapped to those)
nomen_translate = pd.DataFrame([('H0','HS07'), ('H1','HS07'), ('H2','HS07'),
                                ('H3','HS07'), ('H4','HS12'), ('H5','HS12')],
                               columns = ['NomenCode', 'nomen'])

# Add FAO HS nomenclature labels to MFN data
mfn_sub = mfn_sub.merge(nomen_translate, on = 'NomenCode', how = 'left', validate='m:1', indicator = True)
if mfn_sub['_merge'].unique().tolist() != ['both']:
    raise ValueError('Unmatched NomenCodes in raw MFN data.')
else:
    mfn_sub.drop(['_merge'], axis = 1, inplace = True)

# Add ITPD sector codes to MFN data
mfn_sub = mfn_sub.merge(itpd_hs, how ='left', left_on = ['ProductCode', 'nomen'], right_on = ['hs6', 'nomen'])

# Aggregate tariffs by Item Code (simple average over tariff lines within hs6 within itpd Code)
mfn_sub = mfn_sub[['Reporter_ISO_N', 'Year','itpd','TotalNoOfLines', 'Sum_Of_Rates']]
mfn_fao = mfn_sub.groupby(['Reporter_ISO_N', 'Year','itpd']).sum().reset_index()
mfn_fao['mfn_average'] = mfn_fao['Sum_Of_Rates']/mfn_fao['TotalNoOfLines']
mfn_fao.drop(['Sum_Of_Rates'], axis = 1, inplace = True)


# Add DGD iso country codes
mfn_countries = pd.read_csv(mfn_country_codes_local,  encoding='latin')
mfn_fao = mfn_fao.merge(mfn_countries, how = 'left', left_on = 'Reporter_ISO_N',
                        right_on = 'CountryCode', validate = 'm:1')
mfn_fao = mfn_fao[['ISO3', 'Year', 'itpd', 'TotalNoOfLines', 'mfn_average']]

##
# Expand EUN to EU 28
##

# Load data
eun_membs = pd.read_csv(eun_members_local)
# Add recent years
for yr in ['2017','2018','2019']:
    eun_membs[yr] = eun_membs['2016']

# Reshape Long
eun_membs.set_index(['iso3_eun', 'iso3_memb', 'CountryCode_memb', 'Country_Name_memb'], inplace = True)
eun_membs = eun_membs.stack().reset_index()
eun_membs.rename(columns = {'level_4':'Year',0:'eu_member'}, inplace = True)

# Keep only members in each year
eun_membs = eun_membs.loc[eun_membs['eu_member']==1,:]

# keep only desired years
years_str = [str(year) for year in years]
eun_membs = eun_membs.loc[eun_membs['Year'].isin(years_str),:]

# Recast year as numeric
eun_membs['Year'] = eun_membs['Year'].astype(int)
eun_membs.drop(['CountryCode_memb', 'Country_Name_memb'], axis = 1, inplace= True)

# Split data into EU and not EU frames
eu_tariffs = mfn_fao.loc[mfn_fao['ISO3']=='EUN',:]
not_eu_tariffs = mfn_fao.loc[mfn_fao['ISO3']!='EUN',:]

# Drop HRV in 2013, which still has independent tariff obs. in 2013 but was an EU member for more than half the year
not_eu_tariffs = not_eu_tariffs.loc[(mfn_fao['ISO3']!='HRV')&(mfn_fao['Year']!=2013),:]

# Expand EUN to 28 members
eu_tariffs = eu_tariffs.merge(eun_membs, how = 'outer', left_on = ['ISO3','Year'], right_on=['iso3_eun', 'Year'])

# Clean up merge columns
eu_tariffs.drop(['ISO3', 'iso3_eun','eu_member'], axis = 1, inplace = True)
eu_tariffs.rename(columns = {'iso3_memb':'ISO3'}, inplace= True)

# Recombine
mfn_fao = pd.concat([not_eu_tariffs, eu_tariffs], axis = 0)
mfn_fao.sort_values(['ISO3','Year','itpd'], inplace = True)







# -----
# NTM Data
# -----

# Load data
ntm_raw = pd.read_csv(ntm_local)

# # Generate some summary information
# ntm_countries = ntm_raw['reporter'].unique()
# ntm_yearly_count = ntm_raw[['Year','reporter']].groupby('Year').agg('count')
# ntm_nomen_count = ntm_raw[['NomenCode','reporter']].groupby('NomenCode').agg('count')
#
# # Check for countries with multiple NomenCodes (There were none)
# ntm_nomen_by_country = ntm_raw[['reporter','Year','NomenCode','hs6']].groupby(['reporter','Year','NomenCode']).count().reset_index()
# multi_nomen = ntm_nomen_by_country.loc[ntm_nomen_by_country.duplicated(subset=['reporter','Year']),:]

# Drop discriminatory measures.
ntm_data = ntm_raw.loc[ntm_raw['partner']=='WLD',:].copy()

# Pad hs6 codes with leading zeros
ntm_data['hs6'] = ntm_data['hs6'].astype(str)
ntm_data['hs6'] = ntm_data['hs6'].str.zfill(6)

# # Check match on needed HS codes (all but 1 are present at hs6 level)
# ntm_product_comp = CompareIdentifiers(dataframe_a=ntm_data, code_columns_a='hs6',
#                                       dataframe_b=pd.DataFrame(all_hs, columns = ['hs6']), code_columns_b='hs6')

# Cut down to only the needed HS 6 codes
ntm_data = ntm_data.loc[ntm_data['hs6'].isin(all_hs),:]

# Clean up some unneeded columns
ntm_data.drop(['partner','Partner_ISO_N',"Dataset_id", "NTMNomenclature"], axis = 1, inplace = True)


##
# Match to FAO codes
##

# Convert NomenCode to FAO HS revisions
ntm_data = ntm_data.merge(nomen_translate, how = 'left', left_on = ['NomenCode'], right_on = ['NomenCode'],
                          validate = 'm:1')

# Map to itpd codes
ntm_data = ntm_data.merge(itpd_hs, how = 'left', on = ['hs6','nomen'])
ntm_data.drop(['NomenCode','nomen'], axis = 1, inplace = True)


# Convert to panel:
for year in years:
    ntm_data[str(year)] = 0
    ntm_data.loc[(ntm_data['StartDate']<=year) & (ntm_data['EndDate']>=year),str(year)] = 1

ntm_data = ntm_data.melt(id_vars = ['reporter', 'Reporter_ISO_N', 'hs6', 'ntmcode', 'nbr', 'Year', 'ntm_1_digit',
                                    'StartDate', 'EndDate', 'itpd', 'itpd_description', 'hs6_count'])
ntm_data.rename(columns = {'Year':'obs_year','variable':'year'}, inplace = True)

# Drop years without NTMs
ntm_data = ntm_data.loc[ntm_data['value']==1,:]
ntm_data.drop(['value','StartDate','EndDate'], axis = 1, inplace = True)

ntm_data.sort_values(['reporter', 'year', 'itpd', 'hs6', 'ntmcode'], inplace = True)

# Calculate the total number (nbr) of NTMs by 1-digit NTM code for each country, Item, and year.
#   The mean hs6_count gives the number of HS 6-digit codes in each Item (hs6_code does not vary over group)
agg_ntm_data = ntm_data.groupby(['reporter','itpd','year','Reporter_ISO_N','ntm_1_digit']).agg({'hs6_count':'mean',
                                                                                 'nbr':'sum'}).reset_index()

# Create average number of NTMs per underlying hs6 category
agg_ntm_data['avg_ntm'] = agg_ntm_data['nbr']/agg_ntm_data['hs6_count']

# Recast year as int
agg_ntm_data['year'] =agg_ntm_data['year'].astype(int)


##
# Expand EU
##

eu_ntm_data = agg_ntm_data.loc[agg_ntm_data['reporter']=='EUN',:].copy()
not_eu_ntm_data = agg_ntm_data.loc[agg_ntm_data['reporter']!='EUN',:].copy()

# Combine with Tariff data
eu_ntm_data = eu_ntm_data.merge(eun_membs, how = 'left', left_on = ['reporter','year'], right_on=['iso3_eun', 'Year'])

# Clean up merge columns
eu_ntm_data.drop(['reporter', 'iso3_eun','eu_member', 'Year'], axis = 1, inplace = True)
eu_ntm_data.rename(columns = {'iso3_memb':'reporter'}, inplace= True)

# Recombine
agg_ntm_data = pd.concat([not_eu_ntm_data, eu_ntm_data], axis = 0)
agg_ntm_data.sort_values(['reporter', 'year', 'itpd'], inplace = True)

# Drop unwanted second country identifier
agg_ntm_data.drop(['Reporter_ISO_N'], inplace =True, axis = 1)

# Reshape NTM types wide
agg_ntm_data.set_index(['reporter', 'itpd', 'year', 'ntm_1_digit'], inplace = True)
agg_ntm_data = agg_ntm_data.unstack('ntm_1_digit')
agg_ntm_data.reset_index(inplace = True)

# Flatten multi-index column names
agg_ntm_data.columns = ['_'.join(col).strip('_') for col in agg_ntm_data.columns.values]

# Drop hs count and ntm count columns
hs_count_cols = [col for col in agg_ntm_data.columns if col.startswith('hs6_count') or col.startswith('nbr')]
agg_ntm_data.drop(hs_count_cols, axis = 1, inplace = True)

# Fill missing with zero
average_cols = [col for col in agg_ntm_data.columns if col.startswith('avg_ntm')]
for col in average_cols:
    agg_ntm_data[col] = agg_ntm_data[col].fillna(0)

# Create a list of countries with NTM data
ntm_countries = pd.DataFrame(agg_ntm_data['reporter'].unique())
ntm_countries.to_csv("{}//ntm_data_countries.csv".format(save_local), index = False)
# ----
# GDP Data
# ----
# Reformat to panel shape
gdp_data = pd.read_csv(gdp_local)

# Country Names
gdp_names = gdp_data[['Country Name', 'Country Code']].copy()
gdp_names.drop_duplicates(inplace = True)

# Drop some unwanted columns and reshape years to long format
gdp_data = gdp_data.drop(['Country Name', 'Series Code'], axis = 1)
gdp_data.set_index(['Country Code', 'Series Name'], inplace = True)
gdp_data = gdp_data.stack().reset_index()
gdp_data.rename(columns = {'level_2':'year','Country Code':'iso3_d'}, inplace = True)

# Reshape GDP and GDPPC to wide
gdp_data.set_index(['iso3_d', 'year', 'Series Name'], inplace = True)
gdp_data = gdp_data.unstack(2).reset_index()
gdp_data.columns = ['iso3_d','year','GDP_nom','GDPPC_nom']

# Replace missing value indicators with nan
gdp_data = gdp_data.replace('..',np.nan)

# Clean up year variable and cast as int
gdp_data['year'] = gdp_data['year'].str[0:4].astype(int)


# Replace some ISO3 codes (Kosovo)
gdp_data["iso3_d"].replace({'XKX':'KSV'}, inplace = True)
# Drop Channel Islands (Guernsey + Jersey)
gdp_data = gdp_data.loc[gdp_data['iso3_d']!='CHI',:]

# # Check Country code matching
# grav_sub = pd.read_csv(grav_local[2])
# dgd_codes = grav_sub[['iso3_d','country_d']].copy()
# dgd_codes = dgd_codes.drop_duplicates()
# gdp_check = CompareIdentifiers(dataframe_a=gdp_data, code_columns_a='iso3_d',
#                                 dataframe_b=dgd_codes, code_columns_b='iso3_d')

# ----
# Combine Data Sources
# ----

# Rename some variables so they all match
mfn_fao.rename(columns = {'ISO3':'iso3_d', 'Year':'year', 'TotalNoOfLines':'num_mfn_lines'}, inplace = True)
agg_ntm_data.rename(columns = {'reporter':'iso3_d'}, inplace = True)


unilateral_data = pd.merge(mfn_fao, agg_ntm_data, how = 'outer', on = ['iso3_d', 'year', 'itpd'], validate='1:1')
unilateral_data = unilateral_data.merge(gdp_data, how = 'left', on = ['iso3_d', 'year'], validate='m:1')

unilateral_data = unilateral_data.loc[unilateral_data['year'].isin(years),:]


# Fill missing NTM columns with zero
for col in average_cols:
    unilateral_data[col] = unilateral_data[col].fillna(0)

# Rename ITPD sector code column
unilateral_data.rename(columns = {'itpd':'industry_id'}, inplace = True)
unilateral_data.to_csv("{}\\unilateral_data.csv".format(save_local), index = False)












# ----
# Prep ITPD data
# ----

# # Drop unwanted info 
itpd_raw = pd.read_csv(itpd_local)
itpd_data = itpd_raw.loc[itpd_raw['broad_sector']=='Agriculture',:].copy()
itpd_data.drop(['exporter_name','importer_name'], axis = 1, inplace = True)
itpd_data = itpd_data.loc[itpd_data['year'].isin(years),:].copy()

itpd_data.rename(columns={'exporter_iso3':'exporter', 'importer_iso3':'importer', 'trade':'trade_value'}, inplace = True)
itpd_data.drop(['broad_sector', 'flag_mirror', 'flag_zero'], axis = 1, inplace = True)


# ---
# Gravity Variables
# ---

grav_sets = list()
for path in grav_local:
    grav_load = pd.read_csv(path)
    grav_sets.append(grav_load)
grav_data = pd.concat(grav_sets, axis = 0)
grav_data = grav_data[['iso3_o', 'iso3_d', 'year']+grav_vars]
grav_data = grav_data.loc[grav_data['year'] >= min(years), :]
grav_data['ln_distance'] = np.log(grav_data['distance'])
grav_data['colony_ever'] = grav_data[['colony_of_destination_ever', 'colony_of_origin_ever']].max(axis=1)
grav_data['international'] = 0
grav_data.loc[grav_data['iso3_o'] != grav_data['iso3_d'], 'international'] = 1
grav_data = grav_data.drop(['distance', 'colony_of_destination_ever', 'colony_of_origin_ever'], axis=1)
grav_data.rename(columns={'iso3_o': 'exporter', 'iso3_d': 'importer'}, inplace=True)
grav_data['year'] = grav_data['year']

# Add gdp at the bilateral level
gdp_rename = gdp_data[['iso3_d','year','GDP_nom']].copy()
gdp_rename.rename(columns={'iso3_d':'importer', 'GDP_nom':'GDP_nom_full'}, inplace =True)
grav_data = grav_data.merge(gdp_rename, how = 'left', on = ['importer', 'year'])


est_panel = itpd_data.merge(grav_data, how='left', on=['importer', 'exporter', 'year'])

est_panel.to_csv('{}/itpd_bilateral_panel.csv'.format(save_local), index=False)





