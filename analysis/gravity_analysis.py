__Author__ = "Peter Herman"
__Project__ = "A Pragmatic Approach to Estimating Nondiscriminatory Non-tariff Trade Costs"
__Created__ = "August 06, 2021"
__Description__ = '''First stage gravity estimations amd HLY gravity analysis.'''

import pandas as pd
import os as os
import time as time
from datetime import datetime
import gme as gme
import sys

# Add root directory to system path so that local modules can be located and loaded
sys.path.append("D:\\work\\Peter_Herman\\projects\\estimating_nondiscriminatory_trade_costs")
from economic_analysis_tools.data_analysis.TraderRanking import TraderRanking



# Import and configure Stata
import stata_setup
stata_setup.config("C:\\Program Files\\Stata17", "mp")
from pystata import stata

# ------
# Specification
# ------

root_directory = 'D:\work\Peter_Herman\projects\\estimating_nondiscriminatory_trade_costs\\'


# Years to include
years = list(range(2012,2016))

# Share of total trade value to cover with included countries
trade_coverage = 0.95

# Elasticity of substitution (for HLY standard errors delta method)
sigma = 5


sector_list = list(range(1,27))

# Use subsets of the sectors for queuing up multiple python instances to run simultaneously
# sector_list = list(range(1,4))
# sector_list = list(range(4,7))
# sector_list = list(range(7,10))
# sector_list = list(range(10,13))
# sector_list = list(range(13,16))
# sector_list = list(range(16,19))
# sector_list = list(range(19,23))
# sector_list = list(range(23,27))

save_local = "{}analysis/results".format(root_directory)

# Data Paths 
bilateral_data_local = "{}analysis/constructed_data/itpd_bilateral_panel.csv".format(root_directory)
unilateral_data_local = "{}analysis/constructed_data/unilateral_data.csv".format(root_directory)

ntm_countries_local ='{}analysis\constructed_data\\ntm_data_countries.csv'.format(root_directory)


# ---
# Create Save Directories
# ---
for direct in ['DT2_expend', 'DT2_gdp', 'DA2_expend', 'FT2_expend',  'DT3_expend', 'DT2_HLY']:
    directory_path = os.path.join(save_local, direct)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)


# ---
# Load data
# ---

bilat_data = pd.read_csv(bilateral_data_local)
bilat_sectors = bilat_data['industry_id'].unique().tolist()

bilat_data['clusterid'] = bilat_data['exporter']+bilat_data['importer']

unilat_data = pd.read_csv(unilateral_data_local)

timer_history = [('Specification','Sector', 'Minutes', 'Completion date')]

# ---
# Add NTM and Tariff data
# ---
# Set internal measures equal to zero

unilat_measures = unilat_data[['iso3_d','year', 'industry_id', 'mfn_average',
                               'avg_ntm_A', 'avg_ntm_B', 'avg_ntm_C']].copy()
unilat_measures.rename(columns = {'iso3_d':'importer'}, inplace=  True)
bilat_data = bilat_data.merge(unilat_measures, how = 'left', on = ['importer','year', 'industry_id'])

# Set intranational values to zero
for col in ['mfn_average', 'avg_ntm_A', 'avg_ntm_B', 'avg_ntm_C']:
    bilat_data.loc[bilat_data['importer']==bilat_data['exporter'],col] = 0

# Define a composite measure for all technical measures
bilat_data['avg_ntm_tech'] = bilat_data['avg_ntm_A']+bilat_data['avg_ntm_B']+bilat_data['avg_ntm_C']


# ----
# Get Desired Data Subset
# ----
# Select Years
est_data = bilat_data.loc[bilat_data['year'].isin(years), :].copy()

# Compute top sectors by value
top_items = bilat_data[['trade_value','industry_id']].groupby(['industry_id']).sum()
top_items.sort_values('trade_value', inplace = True, ascending=False)

# Create Trade divided by GDP (rescaled to millions of  dollars to match ITPD flows)
est_data['trade_gdp'] = est_data['trade_value'] / (est_data['GDP_nom_full'] / 1000000)

# Create Expenditure and Trade/Expenditure Data
expenditures = est_data.groupby(['importer', 'year', 'industry_id']).agg({'trade_value':'sum'}).reset_index()
expenditures.rename(columns={'trade_value':'expenditure'}, inplace = True)
est_data = est_data.merge(expenditures, on = ['importer', 'year', 'industry_id'], how = 'left',
                                              validate = 'm:1')
est_data['trade_expend'] = est_data['trade_value']/ est_data['expenditure']



no_domestic = est_data.loc[est_data['importer']!=est_data['exporter'],:].copy()


# Pare down countries in sample for restricted sample
bilat_ranking = TraderRanking(est_data, sector_var_name='industry_id')
top_traders = bilat_ranking.cumulative_coverage(cumulative_percentage=trade_coverage,by_year=False, by_sector=True)
trader_num = [(sect,len(top_traders[sect])) for sect in top_traders.keys()]
trader_num_2 = [len(top_traders[sect]) for sect in top_traders.keys() if sect!=15]
ave_trader = sum(trader_num_2)/len(trader_num_2)
all_traders = set()
for sect in top_traders.keys():
    all_traders = all_traders.union(set(top_traders[sect]))


# Create panel of included countries
included_obs = [(exp, imp, yr, sect) for sect in top_traders.keys() for exp in top_traders[sect] for imp in top_traders[sect] for yr in years]
included_panel = pd.DataFrame(included_obs, columns = ['exporter','importer','year','industry_id'])
top_countries_data = included_panel.merge(est_data, how ='left', on = ['exporter', 'importer', 'year', 'industry_id'], validate ='1:1')
top_countries_data = top_countries_data.loc[~top_countries_data['trade_value'].isna(),:]



top_countries_no_domestic = top_countries_data.loc[top_countries_data['importer']!=top_countries_data['exporter'],:].copy()

# Data with top countries and only importers with NTM information
ntm_countries = pd.read_csv(ntm_countries_local)
ntm_countries = ntm_countries['0'].to_list()
top_countries_ntm_data = top_countries_data.loc[top_countries_data['importer'].isin(ntm_countries),:].copy()

# ---
# Create NTM Summary Infromation
# ---
ntm_data = unilat_data.loc[unilat_data['iso3_d'].isin(ntm_countries),:].copy()
ntm_data['avg_ntm_tech'] = ntm_data['avg_ntm_A']+ntm_data['avg_ntm_B']+ntm_data['avg_ntm_C']

ntm_summary_all = ntm_data[['avg_ntm_tech','year']].groupby('year').describe()
ntm_summary_all['type']='all countries'
# Figure out top countries in each sector
top_importers = [(iso, yr, sect) for sect in top_traders.keys() for iso in top_traders[sect] for yr in years]
top_importers = pd.DataFrame(top_importers, columns = ['iso3_d','year','industry_id'])
top_importers['top'] = 1
# Merge dataframe of top countries with the NTM observations
ntm_top_countries = ntm_data.merge(top_importers, how = 'left', on = ['iso3_d','year','industry_id'])
# Keep only observations corresponding to a top country
ntm_top_countries = ntm_top_countries.loc[ntm_top_countries['top']==1,:]
ntm_summary_top = ntm_top_countries[['avg_ntm_tech','year']].groupby('year').describe()
ntm_summary_top['type']='top countries'

ntm_coverage = pd.concat([ntm_summary_top, ntm_summary_all], axis = 0)
ntm_coverage = ntm_coverage.round(2)
ntm_coverage = ntm_coverage.T
ntm_coverage.to_csv("{}\\summary_ntm_coverage.csv".format(save_local))




def statafy_results(gme_model):
    '''
    Convert GME results to dataframes matching the outputs of the Stata estimations.
    :param gme_model: (gme.EstimationModel) An estimated gme model with only sector 'all'
    :return: (pd.DataFrame, pd.DataFrame) A tuple of dataframes. The first is the parameter estimates, the second is the
    fixed effect estimates
    '''
    results = gme_model.results_dict['all']
    #
    # Build FE output to match stat output
    # Collect estimates for fixed effects
    param_ests = pd.DataFrame(results.params)
    param_ests.reset_index(inplace = True)
    # Get importer fixed effects
    imp_fe_ests = param_ests.loc[param_ests['index'].str.startswith('importer_year'),:].copy()
    imp_fe_ests['importer'] = imp_fe_ests['index'].str[17:20]
    imp_fe_ests['year'] = imp_fe_ests['index'].str[20:25]
    imp_fe_ests.rename(columns = {0:'imp_fe'}, inplace= True)
    imp_fe_ests.drop('index', axis = 1, inplace=True)
    # Get importer fixed effects
    exp_fe_ests = param_ests.loc[param_ests['index'].str.startswith('exporter_year'), :].copy()
    exp_fe_ests['exporter'] = exp_fe_ests['index'].str[17:20]
    exp_fe_ests['year'] = exp_fe_ests['index'].str[20:25]
    exp_fe_ests.rename(columns={0: 'exp_fe'}, inplace=True)
    exp_fe_ests.drop('index', axis=1, inplace=True)
    fe_ests = exp_fe_ests.merge(imp_fe_ests, how = 'outer', on = 'year')
    #
    # Build Params output to match Stata output
    parms = pd.concat([results.params, results.bse, results.pvalues], axis = 1)
    parms.reset_index(inplace=True)
    parms.columns = ['parm','estimate', 'stderr', 'p']
    # Add significance stars
    parms['stars'] = ''
    parms.loc[(parms['p']<=0.10)&(parms['p']>0.05), 'stars'] = '*'
    parms.loc[(parms['p']<=0.05)&(parms['p']>0.01), 'stars'] = '**'
    parms.loc[(parms['p'] <= 0.01), 'stars'] = '***'
    
    parms['eq'] = 'trade_gdp'
    
    # Drop fixed effect estimates
    parms = parms.loc[~(parms['parm'].str.startswith('importer_year')|parms['parm'].str.startswith('exporter_year')),:]
    
    return parms, fe_ests


# ---
# (DT2 expend) domestic trade flows, Top countries, 2-way gravity, trade/expenditure GME
# ---
for sector in sector_list:
    print(" **** {} ****".format(sector))
    tic = time.perf_counter()
    two_way_save = "{}//DT2_expend//{}".format(save_local, sector)
    sector_data = top_countries_data.loc[top_countries_data['industry_id'] == sector, :].copy()
    gme_data = gme.EstimationData(sector_data)
    gme_model = gme.EstimationModel(gme_data,
                                    lhs_var='trade_expend',
                                    rhs_var=['international', 'ln_distance', 'contiguity', 'common_language',
                                             'colony_ever', 'agree_pta', 'member_eu_joint', 'member_wto_joint'],
                                    fixed_effects=[['importer','year'],['exporter','year']],
                                    omit_fixed_effect=[['exporter','year']])
    gme_model.estimate()
    parms, fe_ests = statafy_results(gme_model)
    parms.to_stata("{}_param_ests.dta".format(two_way_save))
    fe_ests.to_csv("{}_fe_ests.csv".format(two_way_save))


# ---
# (DT2 gdp) Domestic trade flows, Top countries, 2-way gravity, trade/GDP GME
# ---
for sector in sector_list:
    print(" **** {} ****".format(sector))
    tic = time.perf_counter()
    two_way_save = "{}//DT2_gdp//{}".format(save_local, sector)
    sector_data = top_countries_data.loc[top_countries_data['industry_id'] == sector, :].copy()
    gme_data = gme.EstimationData(sector_data)
    gme_model = gme.EstimationModel(gme_data,
                                    lhs_var='trade_gdp',
                                    rhs_var=['international', 'ln_distance', 'contiguity', 'common_language',
                                             'colony_ever', 'agree_pta', 'member_eu_joint', 'member_wto_joint'],
                                    fixed_effects=[['importer','year'],['exporter','year']],
                                    omit_fixed_effect=[['exporter','year']])
    gme_model.estimate()
    parms, fe_ests = statafy_results(gme_model)
    parms.to_stata("{}_param_ests.dta".format(two_way_save))
    fe_ests.to_csv("{}_fe_ests.csv".format(two_way_save))


# ---
# (DA2 expend) domestic trade, all countries, 2-way gravity, trade/expenditure, GME
# ---
for sector in sector_list:
    print(" **** {} ****".format(sector))
    tic = time.perf_counter()
    two_way_save = "{}//DA2_expend//{}".format(save_local, sector)
    sector_data = est_data.loc[est_data['industry_id'] == sector, :].copy()
    gme_data = gme.EstimationData(sector_data)
    gme_model = gme.EstimationModel(gme_data,
                                    lhs_var='trade_expend',
                                    rhs_var=['international', 'ln_distance', 'contiguity', 'common_language',
                                             'colony_ever', 'agree_pta', 'member_eu_joint', 'member_wto_joint'],
                                    fixed_effects=[['importer','year'],['exporter','year']],
                                    omit_fixed_effect=[['exporter','year']])
    gme_model.estimate()
    parms, fe_ests = statafy_results(gme_model)
    parms.to_stata("{}_param_ests.dta".format(two_way_save))
    fe_ests.to_csv("{}_fe_ests.csv".format(two_way_save))


# ---
# (FT2 expend) Foreign only, top countries, 2-way gravity, trade/expenditure, GME
# ---
for sector in sector_list:
    print(" **** {} ****".format(sector))
    tic = time.perf_counter()
    two_way_save = "{}//FT2_expend//{}".format(save_local, sector)
    sector_data = top_countries_no_domestic.loc[top_countries_no_domestic['industry_id'] == sector, :].copy()
    gme_data = gme.EstimationData(sector_data)
    gme_model = gme.EstimationModel(gme_data,
                                    lhs_var='trade_expend',
                                    rhs_var=['ln_distance', 'contiguity', 'common_language',
                                             'colony_ever', 'agree_pta', 'member_eu_joint', 'member_wto_joint'],
                                    fixed_effects=[['importer','year'],['exporter','year']],
                                    omit_fixed_effect=[['exporter','year']])
    gme_model.estimate()
    try:
        parms, fe_ests = statafy_results(gme_model)
        parms.to_stata("{}_param_ests.dta".format(two_way_save))
        fe_ests.to_csv("{}_fe_ests.csv".format(two_way_save))
    except:
        print('Estimation failed')



# ----
# (DT3 expend) Domestic flows, top countries, three-way fixed effects, trade/expenditure, ppml_panel_sg
# ----
for sector in sector_list:
    print(" **** {} ****".format(sector))
    tic = time.perf_counter()
    three_way_save = "{}//DT3_expend//{}".format(save_local, sector)
    sector_data = top_countries_data.loc[top_countries_data['industry_id'] == sector, :].copy()
    stata.pdataframe_to_data(sector_data, True)
    # Run Stata commands in Python
    stata.run('''
            ppml_panel_sg trade_expend agree_pta member_eu_joint member_wto_joint, ex(exporter) im(importer) y(year) ///
            cluster(clusterid) genS(exp_fe) genM(imp_fe) maxiter(200000)
            parmest, saving("{}", replace) stars(0.1 0.05 0.01)
            keep exporter importer year exp_fe imp_fe
            export delimited using "{}.csv", replace
            '''.format(three_way_save+'_param_ests', three_way_save+'_fe_ests'), echo=True)
    toc = time.perf_counter()
    print("'{}' took {} minutes to complete. {} sectors remaining.".format(sector, round((toc-tic)/60,1),
                                                                     len(sector_list) - sector_list.index(sector) - 1))
    print("{}\n\n".format('*'*100))
    timer_history.append(('three_way', str(sector), str(round((toc-tic)/60,1)), datetime.now().strftime("%d/%m/%Y %H:%M:%S")))



# ----
# (DT2 expend HLY) domestic trade, top countries, 2-way gravity, trade/expenditure, NTM/tariffs, ppml_panel_sg
# ----
for sector in sector_list:
    print(" **** {} ****".format(sector))
    tic = time.perf_counter()
    three_way_save = "{}//DT2_HLY//{}".format(save_local, sector)
    sector_data = top_countries_ntm_data.loc[top_countries_ntm_data['industry_id'] == sector, :].copy()
    stata.pdataframe_to_data(sector_data, True)
    # Run Stata commands in Python
    stata.run('''
            capture log close
            log using "{}.log", replace
            tostring year, gen(year_str)
            gen imp_year = importer+ "_" + year_str
            gen exp_year = exporter+ "_" + year_str
            ppmlhdfe trade_value international ln_distance contiguity common_language colony_ever agree_pta ///
                member_eu_joint member_wto_joint avg_ntm_tech mfn_average, ///
                absorb(imp_year exp_year) vce(robust) maxiter(500000)
            parmest, saving("{}", replace) stars(0.1 0.05 0.01)
            export delimited using "{}.csv", replace
            nlcom -100*(exp(_b[avg_ntm_tech]/{})-1)
            log close
            '''.format(three_way_save+'_log_with_delta_se', three_way_save+'_param_ests', three_way_save+'_fe_ests', sigma), echo=True)
    toc = time.perf_counter()
    print("'{}' took {} minutes to complete. {} sectors remaining.".format(sector, round((toc-tic)/60,1),
                                                                     len(sector_list) - sector_list.index(sector) - 1))
    print("{}\n\n".format('*'*100))
    timer_history.append(('three_way', str(sector), str(round((toc-tic)/60,1)), datetime.now().strftime("%d/%m/%Y %H:%M:%S")))



print("complete!")

