__Author__ = "Peter Herman"
__Project__ = "A Pragmatic Approach to Estimating Nondiscriminatory Non-tariff Trade Costs"
__Created__ = "August 11, 2021"
__Description__ = '''Second Stage AVE cost analysis.'''


import pandas as pd
import numpy as np
import os as os
import statsmodels.api as sm
from economic_analysis_tools.data_analysis.format_regression_table import format_regression_table
import sys
# Load and configure stata for python
import stata_setup
stata_setup.config("C:\\Program Files\\Stata17", "mp")
from pystata import stata # May erroneously show import error in IDE as pystata is loaded from stata and may not be recognized by some IDEs


# ----
# Specifications
# ----

# Set root directory for relative paths
root_directory = 'D:\work\Peter_Herman\projects\\estimating_nondiscriminatory_trade_costs'
sys.path.append(root_directory)

# Save info
save_local = "{}//analysis/results".format(root_directory)

# Set and Omit Sectors
omit_sectors = [15] # 15 has only 6 countries and very few observations for second stage, which prevents it from estimating
sector_list = [str(num) for num in range(1,27) if num not in omit_sectors]

# Years included
year_list = list(range(2010,2016))

# Elasticity of Substitution
sigma = 5
# Level of significance at which to accept AVE estimates as (consider AVE = 0 if pvalue > 1 - est_significance)
est_significance = 0.90

# Specify Estimates from Gravity Stage
DT2_expend_local = "{}/analysis/results/DT2_expend".format(root_directory)
DT2_gdp_local    = "{}/analysis/results/DT2_gdp".format(root_directory)
DA2_expend_local = "{}/analysis/results/DA2_expend".format(root_directory)
FT2_expend_local = "{}/analysis/results/FT2_expend".format(root_directory)
DT3_expend_local = "{}/analysis/results/DT3_expend".format(root_directory)
DT2_expend_HLY_local = "{}/analysis/results/DT2_HLY".format(root_directory)


# Load data
unilateral_data_local = "{}/analysis/constructed_data/unilateral_data.csv".format(root_directory)
bilateral_data_local = "{}/analysis/constructed_data/itpd_bilateral_panel.csv".format(root_directory)
gdp_data_local = "{}\data\source_data\WDI GDP\WDI GDP data 08-05-2021.csv".format(root_directory)
ntm_countries_local = "{}\\analysis/constructed_data/ntm_data_countries.csv".format(root_directory)
sector_labels_local = "{}\\data/source_data/itpd_e_r01/ag_sector_labels.csv".format(root_directory)

# -----
# Load Estimates
# -----


# load bilateral data
bilat_data = pd.read_csv(bilateral_data_local)


##
# Load two-way estimates without NTMs
##

grav_setup = [('DT2_expend', DT2_expend_local),
              ('DT2_gdp', DT2_gdp_local),
              ('DA2_expend', DA2_expend_local),
              ('FT2_expend',FT2_expend_local),
              ('DT3_expend', DT3_expend_local),
              ('DT2_HLY', DT2_expend_HLY_local)]

# Create an accounting of which sector estimates were successfully completed
exist_check = list()
# Dict to store loaded values
loaded_results_dict = dict()

for keyname, grav_dir in grav_setup:
    all_files = os.listdir(grav_dir)
    param_list = list()
    fe_list = list()

    # Load each output file from the gravity estimates
    for file in all_files:
        # Grab sector name
        sector = file.split('_')[0]

        # Load paramater estimates
        if file.endswith('param_ests.dta'):
            loaded = pd.read_stata(os.path.join(grav_dir,file))
            loaded['industry_id'] = sector
            param_list.append(loaded)
            exist_check.append((sector, keyname,'param'))

        # Load Fixed effect estimates
        if file.endswith('fe_ests.csv'):
            loaded = pd.read_csv(os.path.join(grav_dir,file))
            loaded['industry_id'] = sector
            fe_list.append(loaded)
            exist_check.append((sector, keyname,'fe'))

    # Store loaded param values and fixed effect estimates in sub dictionary and then main dictionary
    loaded_dict = dict()
    loaded_dict['param'] = pd.concat(param_list)
    loaded_dict['fe'] = pd.concat(fe_list)
    loaded_results_dict[keyname] = loaded_dict

# Create a matrix of sectors/specifications that were succesfully loaded
existance = pd.DataFrame(exist_check, columns = ['sector', 'specification', 'type'])
existance['exists'] = 1
existance.set_index(['sector', 'specification', 'type'], inplace = True)
existance = existance.unstack('specification')


# Prep Unilateral Data
unilat_data = pd.read_csv(unilateral_data_local)
unilat_data = unilat_data[['iso3_d','year','industry_id', 'GDPPC_nom', 'avg_ntm_A', 'avg_ntm_B', 'avg_ntm_C', 'mfn_average']].copy()
unilat_data['industry_id'] = unilat_data['industry_id'].astype(str)

ntm_countries = pd.read_csv(ntm_countries_local)
ntm_countries = ntm_countries['0'].tolist()
unilat_data = unilat_data.loc[unilat_data['iso3_d'].isin(ntm_countries),:]


# Define the composite measure for all technical measures
unilat_data['avg_ntm_tech'] = unilat_data['avg_ntm_A']+unilat_data['avg_ntm_B']+unilat_data['avg_ntm_C']

# Generate log values
unilat_data['GDPPC_log'] = np.log(unilat_data['GDPPC_nom'])

# Get WTO membership
wto_memb = bilat_data[['importer','year','member_wto_joint']].copy()
wto_memb.rename(columns = {'member_wto_joint':'member_wto', 'importer':'iso3_d'}, inplace = True)
wto_memb = wto_memb.groupby(['iso3_d','year']).max().reset_index()
unilat_data = unilat_data.merge(wto_memb,how = 'left', on = ['iso3_d','year'])

# Get production numbers
prod_data = bilat_data[['exporter','year','industry_id', 'trade_value']].copy()
prod_data = prod_data.groupby(['exporter','year','industry_id']).sum().reset_index()
prod_data.rename(columns = {'trade_value':'production', 'exporter':'iso3_d'}, inplace = True)
prod_data['industry_id'] = prod_data['industry_id'].astype(str)
# Add production to unilateral data
unilat_data = unilat_data.merge(prod_data, how = 'left', on = ['iso3_d','year','industry_id'])
# Fill in cases where countries have no positive exports (domestic or international) with zero production
unilat_data['production'].fillna(0, inplace = True)
# Create log value of production (add 1 dollar to get rid of zeros, which is larger than the min production value)
unilat_data['log_production'] = np.log(unilat_data['production']+(1/1000000))

# ---
# Create IVs for IV version
# ---

##
# Determine 5 closest countries
##
closest_countries = bilat_data.loc[bilat_data['year']==2014,:].copy()
# Drop countries that are problematic neighbors due to data availability
dont_use_neighbors = ['PRK','MAF','NCL','GIB','PYF','VGB']
closest_countries.loc[closest_countries['exporter'].isin(dont_use_neighbors)]
# Drop multiple obs over different sectors and intranational
closest_countries = closest_countries.groupby(['exporter','importer']).max('ln_distance').reset_index()
closest_countries = closest_countries.loc[closest_countries['importer']!=closest_countries['exporter'],
                                          ['importer','exporter', 'ln_distance']]
# Sort by distance and keep 5 shortest by importer
closest_countries.sort_values(['importer','ln_distance'], inplace = True)
closest_countries = closest_countries.groupby('importer').head(5)

##
# Prep and add GDP data
##
# Select necessary data series
gdp_data = pd.read_csv(gdp_data_local)
gdp_data = gdp_data.loc[gdp_data['Series Name']=='GDP (current US$)',:]
gdp_data.set_index(['Country Code'], inplace = True)
gdp_data.drop(['Country Name', 'Series Name','Series Code'], axis = 1, inplace = True)
gdp_data.columns = [int(year[0:4]) for year in gdp_data.columns]
# Use prior year to fill missing observations
for year in gdp_data.columns:
    if year > 2005:
        gdp_data.loc[gdp_data[year]=='..', year] = gdp_data[year-1]
# Reformat long
gdp_data = gdp_data.stack().reset_index()
gdp_data.columns = ['exporter','year', 'exporter_gdp']
gdp_data['year'] = gdp_data['year'].astype(int)
# Convert GDP to numeric
gdp_data.loc[gdp_data['exporter_gdp']=='..','exporter_gdp'] = np.nan
gdp_data['exporter_gdp'] = gdp_data['exporter_gdp'].astype(float)
# Fill remaining missing GDPs (cases missing initial years) with the minimum value for each country
gdp_data["exporter_gdp"] = gdp_data[['exporter',"exporter_gdp"]].groupby("exporter").\
    transform(lambda x: x.fillna(x.min()))
gdp_data = gdp_data.loc[gdp_data['year'].isin(year_list),:]


# Add GDP data to closest countries
closest_countries = closest_countries[['importer','exporter']].merge(gdp_data, on='exporter',
                                                                     how = 'left', validate ='m:m')

# compute total GDP among 5 closest neighbors per importer
total_gdp_neighbors = closest_countries[['importer','year','exporter_gdp']].groupby(['importer','year']).sum()
total_gdp_neighbors.reset_index(inplace = True)
total_gdp_neighbors.rename(columns = {'exporter_gdp':'exporter_gdp_total'}, inplace = True)
closest_countries = closest_countries.merge(total_gdp_neighbors, how = 'left', on = ['importer','year'])

##
# Compute weighted NTM values
##
# Add NTM data
closest_countries = closest_countries.merge(unilat_data[['iso3_d','year','industry_id','avg_ntm_tech','mfn_average']],
                                            how = 'left', left_on = ['exporter','year'], right_on = ['iso3_d','year'])
# Compute gdp weighted ntms (computed as sum[NTM_i*GDP_i/sum(GDP)] := sum(NTM_i*GDP_i)/sum(GDP_i))
closest_countries['closest_ntms_tech'] = closest_countries['exporter_gdp']*\
                                      closest_countries['avg_ntm_tech']/closest_countries['exporter_gdp_total']
# Compute gdp weighted mfn tariffs (computed as sum[MFN_i*GDP_i/sum(GDP)] := sum(MFN_i*GDP_i)/sum(GDP_i))
closest_countries['closest_mfns'] = closest_countries['exporter_gdp']*\
                                      closest_countries['mfn_average']/closest_countries['exporter_gdp_total']
closest_ntms = closest_countries.groupby(['importer','year','industry_id']).agg({'closest_ntms_tech':'sum','closest_mfns':'sum'}).reset_index()


# Rename columns and add to unilat_data
closest_ntms.rename(columns = {'importer':'iso3_d'}, inplace = True)
unilat_data = unilat_data.merge(closest_ntms, how = 'left', on = ['iso3_d', 'year','industry_id'], validate = '1:1')


##
# Compute Total exports and change in imports
##
# Git trade values
trade_values = bilat_data[['exporter','importer','year','industry_id','trade_value']].copy()
# Drop intranational trade
trade_values = trade_values.loc[trade_values['importer']!=trade_values['exporter'],:]
trade_values['industry_id'] = trade_values['industry_id'].astype(str)

# compute total exports
exports = trade_values.groupby(['exporter','year','industry_id']).sum().reset_index()
exports.rename(columns = {'exporter':'iso3_d','trade_value':'total_exports'}, inplace = True)
unilat_data = unilat_data.merge(exports, how = 'left', on = ['iso3_d', 'year','industry_id'], validate = '1:1')

# Compute total imports per importer and industry
imports = trade_values.groupby(['importer','year','industry_id']).sum().reset_index()
# Reshape wide
imports.set_index(['importer','industry_id','year'], inplace = True)
imports = imports.unstack('year')
# Fill missing with zero
imports = imports.fillna(0)

# Calculate change in import between the preceding 2 years
for year in year_list[2:]:
    imports[('import_change',year)]= imports[('trade_value',year-1)] - imports[('trade_value',year-2)]
imports = imports.stack().reset_index()
imports.drop(['trade_value'], axis = 1, inplace=True)
imports.rename(columns = {'importer':'iso3_d'}, inplace=True)
unilat_data = unilat_data.merge(imports, how = 'left', on = ['iso3_d', 'year', 'industry_id'], validate = '1:1')



# ---
# Compute two stage AVEs
# ---


def compute_aves_and_covariates(loaded_results_dict, specification, unilat_data, lower_extreme_drop=None,
                                upper_extreme_drop=None ):
    '''
    Compute AVEs for two-stage approach.
    :param loaded_results_dict: (Dict[SlimResults]) Dictionary of loaded results from the gravity estimations.
    :param specification: (str) The specification to use. From grav_setup variable above.
    :param unilat_data: (DataFrame) Unilateral data to combine with AVE estimates to produce second stage sample.
    :param lower_extreme_drop: (float) Lower Percentile to drop fixed effects from ('imp_01_pct', 'imp_05_pct',
        'imp_10_pct')
    :param upper_extreme_drop: (float) Upper percentile to drop fixed effects from ('imp_90_pct', 'imp_95_pct',
        'imp_99_pct')
    :return: A DataFrame of NTM estimates
    '''

    two_stage_params = loaded_results_dict[specification]['param'].copy()
    two_stage_params.sort_values('parm', inplace =True)

    two_stage_fes = loaded_results_dict[specification]['fe'].copy()
    # Grab importer FEs for AVE calculation
    two_stage_fes = two_stage_fes[['industry_id','year','importer', 'imp_fe']].copy()
    two_stage_fes.drop_duplicates(inplace = True)

    ##
    # Compute AVEs
    ##

    # Collect data
    two_stage_years = two_stage_fes['year'].unique()

    # Determine fixed effect percentiles to drop outlying estimates
    quant_fe = two_stage_fes[['industry_id', 'year','imp_fe']].groupby(['industry_id', 'year']).quantile([0.01, 0.05,0.1, 0.90, 0.95, 0.99])
    quant_fe = quant_fe.unstack(2)
    quant_fe.columns = ['imp_01_pct', 'imp_05_pct', 'imp_10_pct', 'imp_90_pct', 'imp_95_pct','imp_99_pct']
    quant_fe.reset_index(inplace = True)

    two_stage_fes = two_stage_fes.merge(quant_fe, how = 'left', on = ['industry_id','year'], validate ='m:1')

    # Drop outliers
    if lower_extreme_drop:
        two_stage_fes = two_stage_fes.loc[(two_stage_fes['imp_fe'] >= two_stage_fes[lower_extreme_drop]),:]
    if upper_extreme_drop:
        two_stage_fes = two_stage_fes.loc[(two_stage_fes['imp_fe'] <= two_stage_fes[upper_extreme_drop]),:]

    # Determine Max FE for baseline
    max_fe = two_stage_fes.groupby(['industry_id', 'year']).agg({'imp_fe':'max'}).reset_index()
    max_fe.rename(columns = {'imp_fe':'max_imp_fe'}, inplace = True)
    two_stage_fes = two_stage_fes.merge(max_fe, how = 'outer', on = ['industry_id', 'year'], validate ='m:1')



    # Compute AVE:  T = 100 * [exp{(mu - mu_benchmark)/(1-sigma)}-1]
    two_stage_fes['fe_diff'] = two_stage_fes['imp_fe'] - two_stage_fes['max_imp_fe']
    ave_computer = lambda x: 100*(np.exp(x / (1 - sigma)) - 1)

    two_stage_fes['total_ave'] = two_stage_fes['fe_diff'].apply(ave_computer)


    # ----
    # Prep Second Stage covariates
    # ----

    # Grab only wanted ave column, drop duplicates due to it originally being bilateral, and merge with unilateral data
    second_data = two_stage_fes[['importer','year','industry_id','total_ave','max_imp_fe']].copy()
    second_data = second_data.loc[~second_data['total_ave'].isna(),:]
    second_data.drop_duplicates(inplace = True)

    second_data = second_data.merge(unilat_data, how = 'left', left_on = ['importer','year','industry_id'],
                                    right_on = ['iso3_d', 'year', 'industry_id'], validate = '1:1')

    # Add exporter FEs
    exp_fes = loaded_results_dict[specification]['fe'].copy()
    # Grab exporter fixed effects to use as covariate
    exp_fes = exp_fes[['industry_id','year','exporter', 'exp_fe']].copy()
    # Prep exporter fixed effects for merge by renaming country col name, dropping bilateral duplicates, and merging
    exp_fes.rename(columns = {'exporter':'importer'},inplace = True)
    # Drop observation without estimates (often due to missing importer fixed effects)
    exp_fes = exp_fes.loc[~exp_fes['exp_fe'].isna(),:]
    exp_fes.drop_duplicates(inplace = True)
    second_data = second_data.merge(exp_fes, how = 'left', on = ['importer','year', 'industry_id'], validate = '1:1')



    # Create constant
    second_data['constant'] = 1

    # Create sector fixed effects
    sector_dummies = pd.get_dummies(second_data['industry_id'], prefix='itpd_fe')
    second_data = pd.concat([second_data, sector_dummies], axis=1)

    # Create year fixed effects
    year_dummies = pd.get_dummies(second_data['year'], prefix='year_fe')
    second_data = pd.concat([second_data, year_dummies], axis=1)

    return second_data

# ---
# Create Second-stage data
# ---
# Data with only top countries + Domestic + two-way gravity
DT2_expend_data = compute_aves_and_covariates(loaded_results_dict, 'DT2_expend', unilat_data)
DT2_expend_sector_dummies = [col for col in DT2_expend_data.columns if col.startswith('itpd_fe')]
year_dummy_names = [col for col in DT2_expend_data.columns if col.startswith('year_fe')]

# Data with domestic flows + Top countries + two way + Trade/GDP
DT2_gdp_data = compute_aves_and_covariates(loaded_results_dict, 'DT2_gdp', unilat_data)
DT2_gdp_sector_dummies = [col for col in DT2_expend_data.columns if col.startswith('itpd_fe')]

# Data with all countries + domestic + two-way
DA2_expend_data = compute_aves_and_covariates(loaded_results_dict, 'DA2_expend', unilat_data)
DA2_expend_sector_dummies = [col for col in DA2_expend_data.columns if col.startswith('itpd_fe')]

# Data with top countries + no domestic + two-way
FT2_expend_data = compute_aves_and_covariates(loaded_results_dict, 'FT2_expend', unilat_data)
FT2_expend_sector_dummies = [col for col in FT2_expend_data.columns if col.startswith('itpd_fe')]

# Data with Domestic + top countries + Three-way gravity
DT3_expend_data = compute_aves_and_covariates(loaded_results_dict, 'DT3_expend', unilat_data,
                                              lower_extreme_drop='imp_05_pct', upper_extreme_drop='imp_95_pct')
DT3_expend_data = DT3_expend_data.loc[DT3_expend_data['industry_id']!='18',:]
DT3_expend_sector_dummies = [col for col in DT3_expend_data.columns if col.startswith('itpd_fe')]



# --
# Run Regressions
# --

# Regression Naming conventions
# (DT2): Domestic trade included, top countries only, 2-way gravity
# (DA2): Domestic trade included, all countries, 2-way gravity
# (FT2): Foreign trade only, top countries, 2-way gravity
# (DT3): Domestic trade included, top countries, three-way

regressions = dict()

##
# DT2: Regressions using only the top countries
##
# Define specifications
DT2_expend_specs = {'(DT2) ntms': ['avg_ntm_tech'] + DT2_expend_sector_dummies + year_dummy_names,
                       '(DT2) controls': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production'] + DT2_expend_sector_dummies + year_dummy_names,
                       '(DT2) max_fe': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production','max_imp_fe'] + DT2_expend_sector_dummies + year_dummy_names,
                       '(DT2) exporter_fe': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'exp_fe','max_imp_fe'] + DT2_expend_sector_dummies + year_dummy_names}
for name, spec in DT2_expend_specs.items():
    model = sm.OLS(endog = DT2_expend_data['total_ave'], exog= DT2_expend_data[spec],
                   missing='drop').fit(cov_type = 'HC1')
    regressions[name] = model

##
# DA2: Regressions that add all countries and have domestic flows
##
DA2_expend_specs = {'(DA2) ntms': ['avg_ntm_tech'] + DA2_expend_sector_dummies,
                       '(DA2) no_max_fe': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production',] + DA2_expend_sector_dummies + year_dummy_names,
                       '(DA2) max_fe': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production','max_imp_fe'] + DA2_expend_sector_dummies + year_dummy_names}

for name, spec in DA2_expend_specs.items():
    model = sm.OLS(endog = DA2_expend_data['total_ave'], exog= DA2_expend_data[spec],
                   missing='drop').fit(cov_type = 'HC1')
    regressions[name] = model

##
# DT3: Domestic, top, and three-way gravity
##
DT3_expend_specs = {'(DT3) production': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production', 'max_imp_fe'] + DT3_expend_sector_dummies + year_dummy_names}

for name, spec in DT3_expend_specs.items():
    model = sm.OLS(endog = DT3_expend_data['total_ave'], exog= DT3_expend_data[spec],
                   missing='drop').fit(cov_type = 'HC1')
    regressions[name] = model
    print(model.summary())


##
# FT2: Regressions with top countries but exclude domestic flows
##
FT2_expend_specs = {'(FT2) production': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production', 'max_imp_fe'] + FT2_expend_sector_dummies + year_dummy_names}

for name, spec in FT2_expend_specs.items():
    model = sm.OLS(endog = FT2_expend_data['total_ave'], exog= FT2_expend_data[spec],
                   missing='drop').fit(cov_type = 'HC1')
    regressions[name] = model



##
# DT2: Trade/GDP
##
DT2_gdp_specs = {'(DT2) trade/gdp': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production', 'max_imp_fe'] + DT2_gdp_sector_dummies + year_dummy_names,
                 '(DT2) trade/gdp exp_fe': ['avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'exp_fe', 'max_imp_fe'] + DT2_gdp_sector_dummies + year_dummy_names}
for name, spec in DT2_gdp_specs.items():
    model = sm.OLS(endog = DT2_gdp_data['total_ave'], exog= DT2_gdp_data[spec],
                   missing='drop').fit(cov_type = 'HC1')
    regressions[name] = model


# ---
# IV estimation
# ---
stata.pdataframe_to_data(DT2_expend_data, True)
# Run Stata commands in Python
iv_save_name = "{}/iv_ests_ntm_mfn".format(save_local)
stata.run('''
        * log close
        log using "{}.log", replace
        // IV for NTMs only
        ivreghdfe total_ave mfn_average member_wto GDPPC_log production max_imp_fe itpd_fe* year_fe* (avg_ntm_tech = closest_ntms_tech total_exports import_change), robust first
        estimates store ntms_only

        // IV for NTMs and MFN
        ivreghdfe total_ave member_wto GDPPC_log production max_imp_fe itpd_fe* year_fe* (avg_ntm_tech mfn_average = closest_ntms_tech closest_mfns total_exports import_change), robust first
        estimates store ntms_mfn

        esttab ntms_only ntms_mfn using "{}.csv", replace star(* 0.1 ** 0.05 *** 0.01) mtitles(ntm_iv ntm_and_mfn_iv) se scalars(F r2_a)
        log close
        '''.format(iv_save_name, iv_save_name), echo=True)

# ---
# Combine all results
# ---



# Create table of preffered specifications
prefered_specs = ['(DT2) ntms', '(DT2) controls', '(DT2) max_fe', '(DT2) exporter_fe', '(DT2) trade/gdp','(DT2) trade/gdp exp_fe', '(DA2) max_fe', '(FT2) production', '(DT3) production']
prefered_regressions = dict()
for key in prefered_specs:
    prefered_regressions[key] = regressions[key]
prefered_results = format_regression_table(prefered_regressions, omit_fe_prefix=['itpd_fe', 'year_fe'], r_squared=True)

# Create new row for adjusted R^2 and fill them
prefered_results.loc['b_adj_R^2', 'Variable'] = 'adj_R^2'
for key in prefered_specs:
    ols_model = regressions[key]
    prefered_results.loc['b_adj_R^2', key] = round(ols_model.rsquared_adj, 3)

# Rename variables
var_rename = [('avg_ntm_tech','NTM'), ('mfn_average','MFN tariff'), ('GDPPC_log','ln(GDPPC)'), ('production','Production'),
 ('member_wto','WTO'), ('exp_fe',"$\\hat{\\mu}_{its}$"), ('max_imp_fe','$\\hat{\\nu}^*_{jts}$'), ('adj_R^2', 'Adj. $R^2$')]
for new_name in var_rename:
    prefered_results = prefered_results.replace(new_name[0], new_name[1])

prefered_results.to_csv("{}//preferred_regressions.csv".format(save_local))
prefered_results.to_csv("{}//preferred_regressions.tex".format(save_local), index = False,
                        sep = "&", line_terminator='\\\\\n')

# Create table of DT2 results
DT2_table = prefered_results[['Variable', '(DT2) ntms', '(DT2) controls', '(DT2) max_fe', '(DT2) exporter_fe', '(DT2) trade/gdp','(DT2) trade/gdp exp_fe']]
DT2_table.to_csv("{}//DT2_table.tex".format(save_local), index = False,
                        sep = "&", line_terminator='\\\\\n')

# Create table of multi-gravity model results
multi_grav_table = prefered_results[['Variable', '(DT2) max_fe', '(DA2) max_fe', '(FT2) production', '(DT3) production']]
multi_grav_table.to_csv("{}//multi_grav_table.tex".format(save_local), index = False,
                        sep = "&", line_terminator='\\\\\n')



# ---
# Compute total NTM costs
# ---
total_ntms = DT2_expend_data[['importer', 'year', 'industry_id', 'total_ave', 'avg_ntm_tech']].copy()
total_ntms = total_ntms.loc[total_ntms['importer'].isin(ntm_countries),:]

ntm_est = regressions['(DT2) max_fe'].params['avg_ntm_tech']
total_ntms['total_ntm'] = total_ntms['avg_ntm_tech']*ntm_est
total_ntms['ntm_share_of_total'] = total_ntms['total_ntm']/total_ntms['total_ave']
total_ntm_all_summary = total_ntms[['total_ntm']].describe()
industry_ntm_summary = total_ntms.groupby(['industry_id'])['total_ntm'].describe().sort_values('mean')
total_ntm_all_summary = pd.concat([total_ntm_all_summary.T, industry_ntm_summary], axis = 0)

country_ntms = total_ntms.groupby(['importer'])['total_ntm'].describe().sort_values('mean').reset_index()
total_ntm_all_summary.to_csv("{}//total_ntm_summary_by_industry.csv".format(save_local))

# Format table for paper and add industry labels
sect_labels = pd.read_csv(sector_labels_local,dtype=str)
industry_ntm_table = total_ntm_all_summary.reset_index()
industry_ntm_table['index'].replace({'total_ntm':'0'}, inplace = True)
industry_ntm_table = industry_ntm_table.merge(sect_labels, left_on='index', right_on ='industry_id', how = 'left')
industry_ntm_table[' industry_name'].fillna('All sectors', inplace = True)
industry_ntm_table['index'] = industry_ntm_table['index'].astype(int)
industry_ntm_table.sort_values('index', inplace = True)
industry_ntm_table = industry_ntm_table[[' industry_name', 'mean', 'std', 'min', '50%', 'max']]
industry_ntm_table = industry_ntm_table.round(2)
industry_ntm_table.to_csv("{}//total_ntm_summary_by_industry.tex".format(save_local), sep='&',
                          line_terminator='\\\\\n', index = False)








# ----
# Compare to HLY
# ----

##
# Two Stage
##
# Define specifications for second stage for both two-way and three-way versions
two_way_comp_spec = ['constant', 'avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production',
                     'max_imp_fe'] + year_dummy_names
three_way_comp_spec = ['constant', 'avg_ntm_tech', 'mfn_average', 'member_wto', 'GDPPC_log', 'production',
                     'max_imp_fe'] + DT3_expend_sector_dummies + year_dummy_names
# Create dicts to store estimates
two_way_sector_results = dict()
three_way_sector_results = dict()

# Run regressions for each sector individually
for sect in sector_list:
    # Grab data for specific sector
    two_way_sector_subset = DT2_expend_data.loc[DT2_expend_data['industry_id'] == sect, :].copy()
    three_way_sector_subset = DT3_expend_data.loc[DT3_expend_data['industry_id'] == sect, :].copy()

    # Define dict to store two way results
    two_way_results = dict()
    # Attempt to run two-way regression
    try:
        two_way_model = sm.OLS(endog=two_way_sector_subset['total_ave'], exog=two_way_sector_subset[two_way_comp_spec],
                       missing='drop').fit(cov_type='HC1')
        two_way_sector_results[sect] = two_way_model
    except:
        print("Sector {} could not estimate".format(sect))


    three_way_results = dict()

    try:
        three_way_model = sm.OLS(endog = three_way_sector_subset['total_ave'], exog= three_way_sector_subset[three_way_comp_spec],
                       missing='drop').fit(cov_type ='HC1')
        three_way_sector_results[sect] = three_way_model
    except:
        print("Sector {} could not estimate".format(sect))


# Coimpile Estimates into a table format with multiple sectors
def format_sector_estimates(results_dict, spec_name):
    '''
    Format second stage regression results for each sector into a single columns with estimates and standard errors.
    :param results_dict: (dict) A dictionary of OLS regression results in which the keys are industry_id numbers.
    :param spec_name: (str) Name to use for column of estimates.
    :return: (pd.DataFrame) A table of formatted estimates
    '''
    sector_ests = list()
    for sector, model in results_dict.items():
        # Grab desired values
        ests = pd.concat([model.params, model.bse, model.pvalues], axis = 1)
        ests.columns = ['ave_estimate', 'std_err', 'pvalue']
        # Select only ntm estimates
        ests = ests.loc[ests.index == 'avg_ntm_tech',:]
        # Add sector identifier columns
        ests['industry_id'] = sector
        # Create column of significance stars
        ests['stars'] = ''
        ests.loc[ests['pvalue'] <= 0.1, 'stars'] = ests['stars'] + '*'
        ests.loc[ests['pvalue'] <= 0.05, 'stars'] = ests['stars'] + '*'
        ests.loc[ests['pvalue'] <= 0.01, 'stars'] = ests['stars'] + '*'
        # Format strings of estiamtes and standard errors
        ests['ave'] = ests['ave_estimate'].round(3).astype(str)+ests['stars']
        ests['std_err'] = '(' + ests['std_err'].round(3).astype(str) + ')'
        sector_ests.append(ests)
    # Combine sectors
    all_sector_ests = pd.concat(sector_ests, axis = 0)
    # Reformat long
    all_sector_ests.set_index('industry_id', inplace = True)
    all_sector_ests = all_sector_ests[['ave','std_err']].stack()
    all_sector_ests = all_sector_ests.reset_index()
    # Add column name
    all_sector_ests.rename(columns={0:spec_name}, inplace = True)
    return all_sector_ests

two_way_sector_ests = format_sector_estimates(results_dict = two_way_sector_results, spec_name = 'two_way_two_stage')
three_way_sector_ests = format_sector_estimates(results_dict = three_way_sector_results, spec_name = 'three_way_two_stage')

##
# HLY
##


hly_ests = loaded_results_dict['DT2_HLY']['param']
hly_ests.sort_values('parm', inplace = True)
# Define funciton to retrieve ave and delta method standard errors from stata log file.
def log_delta_se_scrape(log_directory):
    '''
    Function to parse delta method standard errors from stata log file.
    :param log_directory: (str) Directory containing the log files.
    :return: (pd.DataFrame) A Dataframe of estimates.
    '''
    log_files = [file for file in os.listdir(log_directory) if file.endswith('delta_se.log')]
    delta_std_errs_list = list()
    for log in log_files:
        with open(os.path.join(log_directory, log)) as log_file:
            sector_id = log.split('_')[0]
            # read file
            log_contents = log_file.readlines()
            # Find line with delta estimates
            std_err = [line for line in log_contents if ('_nl_1' in line) and ('_b[avg_ntm_tech]') not in line]
            # If found, process the line into a usable format
            if std_err != []:
                # Convert a list of what should be only one line into that line as a string.
                std_err = std_err[0]
                # Replace the Stata var name with the sector id
                std_err = std_err.replace('_nl_1',sector_id)
                # Remove a special visual sepperator and split along white spaces.
                std_err = std_err.replace('|',' ')
                std_err = std_err.split()
                delta_std_errs_list.append(std_err)
            else:
                # If there is no estimate in the log, record 'no estimate' for the sector.
                delta_std_errs_list.append([sector_id]+[np.nan]*6)
    # Convert to dataframe and label columns
    delta_df = pd.DataFrame(delta_std_errs_list)
    delta_df.columns = ['industry_id', 'ave_estimate', 'std_err', 'delta_z', 'delta_p', 'delta_95_lower', 'delta_95_upper']
    # Replace missing/omitted estimates with np.nan
    delta_df = delta_df.replace('.',np.nan)
    delta_df = delta_df.replace('(omitted)',np.nan)

    for col in ['ave_estimate', 'std_err', 'delta_z', 'delta_p', 'delta_95_lower', 'delta_95_upper']:
        delta_df[col] = delta_df[col].astype(float)
    delta_df['stars'] = ''
    delta_df.loc[delta_df['delta_p']<=0.1,'stars'] = delta_df['stars']+'*'
    delta_df.loc[delta_df['delta_p'] <= 0.05, 'stars'] = delta_df['stars'] + '*'
    delta_df.loc[delta_df['delta_p'] <= 0.01, 'stars'] = delta_df['stars'] + '*'
    delta_df['ave'] = delta_df['ave_estimate'].round(3).astype(str)+delta_df['stars']
    delta_df['std_err'] = '(' + delta_df['std_err'].round(3).astype(str) + ')'
    return delta_df

# Retrieve standard error estimates from the log files
two_way_delta_se = log_delta_se_scrape(DT2_expend_HLY_local)

# Reformat values into table format
two_way_delta_se.set_index('industry_id', inplace = True)
two_way_delta_se = two_way_delta_se[['ave','std_err']].stack()
two_way_delta_se = two_way_delta_se.reset_index()
two_way_delta_se.rename(columns={0:'two_way_HLY'}, inplace = True)


##
# Combine
##
# Merge 2-stage and HLY
comparrison_ests = two_way_sector_ests.merge(two_way_delta_se, on = ['industry_id','level_1'], how = 'outer')
# Reshape Standard errors to their own columns
comparrison_ests.set_index(['industry_id','level_1'], inplace = True)
comparrison_ests = comparrison_ests.unstack('level_1')
comparrison_ests.reset_index(inplace = True)
comparrison_ests.columns = [col[0]+col[1] for col in comparrison_ests.columns]
# Add labels to sectors
comparrison_ests = sect_labels.merge(comparrison_ests, on = 'industry_id')


# Save csv of values
comparrison_ests.to_csv('{}//2_stage_HLY_comparison.csv'.format(save_local), index = False)
comparrison_ests['industry_id'] = comparrison_ests['industry_id'].astype(int)
comparrison_ests.sort_values(['industry_id'], inplace = True)
comparrison_ests.fillna('-', inplace = True)
comparrison_ests = comparrison_ests.replace('nan','-')
comparrison_ests = comparrison_ests.replace('(nan)','-')
comparrison_ests.drop(['industry_id'], axis = 1, inplace = True)
comparrison_ests.to_csv('{}//2_stage_HLY_comparison.tex'.format(save_local), index = False, sep = '&',
                        line_terminator='\\\\\n')


# ---
# Examine differences in AVE estimates across gravity models
# ---
ave_descriptions = list()
for name, data in [('DT2', DT2_expend_data), ('DA2', DA2_expend_data), ('FT2', FT2_expend_data),
                   ('DT3', DT3_expend_data)]:
    summary = pd.DataFrame(data[['total_ave']].describe())
    summary.columns = [name+'_total_ave']
    ave_descriptions.append(summary)
ave_descriptions = pd.concat(ave_descriptions,axis = 1)
ave_descriptions = ave_descriptions.round(2)
ave_descriptions.to_csv("{}\\ave_summary_stats.tex".format(save_local), sep = '&', line_terminator='\\\\\n')
comp = pd.merge(DT2_expend_data, FT2_expend_data, how ='outer', on = ['importer', 'year', 'industry_id'],
                indicator=True)

# ----
# Data Summary Information
# ----

# Correlation between production and exporter fe
DT2_expend_data[['production', 'exp_fe']].corr()
DT2_expend_data[['year', 'avg_ntm_tech']].groupby('year').describe().T
DA2_expend_data[['year', 'avg_ntm_tech']].groupby('year').describe().T
