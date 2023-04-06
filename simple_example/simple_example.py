__Author__ = "Peter Herman"
__Created__ = "April 06, 2023"
__Description__ = '''A simple example demonstrating the two-stage methodology. The analysis runs the first and 
second stage for wheat trade and a small sample of major importers/exporters. The analysis does not necessarily follow 
all recommendations from the paper for the sake of simplicity. All required data files are available at 
https://doi.org/10.7910/DVN/AGG5T9.'''



import pandas as pd
import gme as gme
import numpy as np
import statsmodels.api as sm



# ------
# Step 1: Load and prep data (Data files available from  https://doi.org/10.7910/DVN/AGG5T9)
# ------

# ITPD bilateral panel
bilat_data = pd.read_csv("analysis/constructed_data/itpd_bilateral_panel.csv")

# Country-level data
unilat_data = pd.read_csv("analysis/constructed_data/unilateral_data.csv")


# Simplify to a subsample of countries
top_countries = ['USA', 'JPN', 'DEU', 'FRA', 'CHN', 'IND', 'ITA', 'ESP', 'NLD', 'KOR', 'MEX', 'CHE', 'BEL', 'HKG',
                 'GBR', 'POL', 'AUT', 'NGA', 'CAN', 'DNK', 'CZE', 'ISR', 'CHL', 'RUS', 'FIN', 'SGP', 'IDN', 'ROU',
                 'BRA', 'PRT', 'IRL', 'MYS', 'AUS', 'TUR', 'THA', 'VNM', 'PER', 'GRC', 'LUX', 'HUN', 'SWE', 'PAK',
                 'SVK', 'SAU', 'UKR', 'PHL', 'ZAF', 'NOR', 'BGR', 'IRN', 'ECU', 'KAZ', 'ARG', 'SVN', 'LKA', 'CRI',
                 'EGY', 'LTU', 'BGD']
bilat_data = bilat_data.loc[bilat_data['exporter'].isin(top_countries) & bilat_data['importer'].isin(top_countries),:]

# Divide Trade flows by importer GDP (divide GDP by million to get $M)
bilat_data['trade_gdp'] = bilat_data['trade_value'] / (bilat_data['GDP_nom_full'] / 1000000)

# Select just wheat trade
wheat_data = bilat_data.loc[bilat_data['industry_descr']=='Wheat',:]

# ----
# Step 1: Estimate 1st stage (using GME package)
# ----

# Create a GME data model
gme_data = gme.EstimationData(wheat_data, imp_var_name='importer', exp_var_name='exporter', year_var_name='year',
                              trade_var_name='trade_value', sector_var_name='industry_descr')

# Define GME estimation specification
est_model = gme.EstimationModel(gme_data,
                                lhs_var='trade_gdp',
                                rhs_var=['contiguity', 'common_language', 'agree_fta', 'member_eu_joint',
                                         'member_wto_joint', 'ln_distance', 'colony_ever', 'international'],
                                sector_by_sector=True,
                                fixed_effects=[['exporter', 'year'], ['importer', 'year']],
                                omit_fixed_effect=[['exporter','year']])

results = est_model.estimate()
# Print estimation results
results['Wheat'].summary()

'''
"""
                 Generalized Linear Model Regression Results                  
==============================================================================
Dep. Variable:              trade_gdp   No. Iterations:                     14
Model:                            GLM   Df Residuals:                     9913
Model Family:                 Poisson   Df Model:                          695
Link Function:                    log   Scale:                          1.0000
Method:                          IRLS   Log-Likelihood:                -8.1723
Covariance Type:                  HC1   Deviance:                      0.60534
No. Observations:               10609   Pearson chi2:                     3.34
============================================================================================
                               coef    std err          t      P>|t|      [0.025      0.975]
--------------------------------------------------------------------------------------------
contiguity                   1.0272      0.116      8.870      0.000       0.800       1.254
common_language              0.2496      0.108      2.303      0.021       0.037       0.462
agree_fta                    0.9625      0.102      9.407      0.000       0.762       1.163
member_eu_joint              1.2012      0.205      5.858      0.000       0.799       1.603
member_wto_joint             1.6285      0.347      4.694      0.000       0.948       2.309
ln_distance                 -1.2282      0.077    -15.930      0.000      -1.379      -1.077
colony_ever                 -0.1582      0.292     -0.542      0.588      -0.731       0.414
international               -5.5288      0.201    -27.469      0.000      -5.923      -5.134
exporter_year_fe_ARG2010     3.8268      1.056      3.623      0.000       1.756       5.897
exporter_year_fe_ARG2011     5.5594      0.920      6.040      0.000       3.755       7.364
exporter_year_fe_ARG2012     5.8551      0.955      6.131      0.000       3.983       7.727
[Other fixed effects omitted for brevity]
'''



# ----
# Step 2: Generate AVEs from importer fixed effects
# ----

# Extract fixed effect estimates from gravity results
parameter_ests = results['Wheat'].params
imp_fe = parameter_ests.loc[parameter_ests.index.str.startswith(('importer_year_fe'))]

# Convert series to data frame and add some identifier columns back in
imp_fe = pd.DataFrame(imp_fe, columns=['imp_fe'])
imp_fe['iso3_d'] = imp_fe.index.str[17:20]
imp_fe['year'] = imp_fe.index.str[20:24]
imp_fe['year'] = imp_fe['year'].astype(int)
imp_fe['industry_descr'] = 'Wheat'
imp_fe['industry_id'] = 1

# Determine baseline/max FE values
max_fe = imp_fe.groupby(['industry_descr', 'year']).agg({'imp_fe':'max'}).reset_index()
max_fe.rename(columns = {'imp_fe':'max_imp_fe'}, inplace = True)
imp_fe = imp_fe.merge(max_fe, how = 'outer', on = ['industry_descr', 'year'], validate ='m:1')

# Compute AVE:  T = 100 * [exp{(mu - mu_benchmark)/(1-sigma)}-1]
imp_fe['fe_diff'] = imp_fe['imp_fe'] - imp_fe['max_imp_fe']
sigma = 7
ave_function = lambda x: 100*(np.exp(x / (1 - sigma)) - 1)
imp_fe['total_ave'] = imp_fe['fe_diff'].apply(ave_function)

# Print stats about total estimated AVEs
imp_fe['total_ave'].describe()


# -----
# Step 4: Prepare Second stage, unilateral data
# -----

ave_data = imp_fe.merge(unilat_data, how = 'left', on=['iso3_d', 'year', 'industry_id'], validate = '1:1')

# Generate year fixed effects
year_fe = pd.get_dummies(ave_data['year'], prefix='year_fe')

ave_data = pd.concat([ave_data, year_fe], axis = 1)

# Define set of controls/covariates (excluding 2010 year fixed effect)
control_vars = ['mfn_average', 'avg_ntm_A', 'GDP_nom', 'GDPPC_nom', 'max_imp_fe', 'year_fe_2011', 'year_fe_2012',
                'year_fe_2013', 'year_fe_2014', 'year_fe_2015',]

second_stage_ests = sm.OLS(endog = ave_data['total_ave'], exog= ave_data[control_vars],
                                  missing='drop').fit(cov_type = 'HC1')
second_stage_ests.summary()

"""
                            OLS Regression Results                            
==============================================================================
Dep. Variable:              total_ave   R-squared:                       0.296
Model:                            OLS   Adj. R-squared:                  0.269
Method:                 Least Squares   F-statistic:                       nan
Date:                Thu, 06 Apr 2023   Prob (F-statistic):                nan
Time:                        15:04:47   Log-Likelihood:                -1451.2
No. Observations:                 245   AIC:                             2922.
Df Residuals:                     235   BIC:                             2957.
Df Model:                           9                                         
Covariance Type:                  HC1                                         
================================================================================
                   coef    std err          z      P>|z|      [0.025      0.975]
--------------------------------------------------------------------------------
mfn_average     -0.7468      0.299     -2.494      0.013      -1.334      -0.160
avg_ntm_A        0.8941      0.226      3.961      0.000       0.452       1.337
GDP_nom       1.045e-11   2.15e-12      4.851      0.000    6.23e-12    1.47e-11
GDPPC_nom        0.0015      0.000      4.788      0.000       0.001       0.002
max_imp_fe      47.1809      8.097      5.827      0.000      31.312      63.050
year_fe_2011    -9.7437     15.722     -0.620      0.535     -40.557      21.070
year_fe_2012    -8.1377     25.035     -0.325      0.745     -57.206      40.931
year_fe_2013   -15.5335     23.802     -0.653      0.514     -62.184      31.117
year_fe_2014    -5.2508     20.493     -0.256      0.798     -45.416      34.914
year_fe_2015    30.5118     15.920      1.917      0.055      -0.691      61.714
==============================================================================
Omnibus:                       80.543   Durbin-Watson:                   1.567
Prob(Omnibus):                  0.000   Jarque-Bera (JB):              254.118
Skew:                           1.400   Prob(JB):                     6.59e-56
Kurtosis:                       7.130   Cond. No.                     1.81e+13
==============================================================================
Notes:
[1] Standard Errors are heteroscedasticity robust (HC1)
[2] The condition number is large, 1.81e+13. This might indicate that there are
strong multicollinearity or other numerical problems.
"""


