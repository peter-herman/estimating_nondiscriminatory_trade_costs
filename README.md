# README 
This repository contains the instructions, code files, and many of the data inputs needed to replicate the analysis presented in "A Pragmatic Approach top Estimating Nondiscriminatroy Non-tariff Trade Costs" by Peter Herman.

### Citation 
Herman, P. R. (2022). A pragmatic approach to estimating nondiscriminatory non-tariff trade costs. Review of International Economics, 1--30. https://doi.org/10.1111/roie.12604


## Data inputs

The analysis is based on multiple source data files, which are listed below. Many of the smaller data sources are included in the git repository at the location specified but the largest of the files must be downloaded from the original data sources. 

The output data from the data construction step are available in a repository at the Harvard Dataverse at https://doi.org/10.7910/DVN/AGG5T9. Steps 2 and 3, which conduct the analysis from sections 3 and 4 of the paper, can be completed using the data available from Dataverse.

Trade Information 
* **ITPD bilateral trade data**: "data\source_data\itpd_e_r01\ITPD_E_R01.csv" (download from origin)
* **ITPD sector concordance (ITPD to FAO FCL)**: "data\source_data\itpd_e_r01\ITPD Agriculture FCL concordance.txt" (included)
* **FAO to HS concordance**: "data\source_data\\FAOSTAT_expanded_item_codes_8-19-2021.csv" (included)
* **Sector labels**: "data/source_data/itpd_e_r01/ag_sector_labels.csv" (included)

MFN Tariffs
* **MFN Data (zipped folder of data files)**: "data\source_data\TRAINS Bulk MFN 2005-2019" (download from origin)
* **MFN tariff country codes**: "data\source_data\TRAINS Bulk MFN 2005-2019\\00_TRAINS_country_codes.csv" (included)

NTM Information
* **NTM Data**: 'data\source_data/ntm_hs6_2010_2018v.12/NTM_hs6_2010_2018 v.12.csv' (download from origin)

Country Information
* **EUN membership**: "data\source_data\EUN_member_codes.csv" (included)
* **DGD gravity files**: (download from origin)
    * "data/source_data/dgd_2.1/release_2.1_2005_2009.csv"
    * "data/source_data/dgd_2.1/release_2.1_2010_2014.csv"
    * "data/source_data/dgd_2.1/release_2.1_2015_2019.csv"
* **GDP Data (more recent than gravity data)**: "data\source_data\WDI GDP\WDI GDP data 08-05-2021.csv" (included)


---
## Replication instructions
The analysis can be conducted by completing running the following three scripts in sequence. Each step is detailed below.

1. **Create data**: analysis/data_construct.py
2. **Conduct first stage gravity analysis:** analysis/gravity_analysis.py
3. **Conduct second stage analysis**: analysis/ave_analysis.py 

### Step 1: Create data
* **Script**: analysis/data_construct.py

This script loads the source data, cleans, concords, and merges the data into two main samples: on bilateral for gravity use and one unilateral for use in both stages. 

#### Data Inputs

Trade Information
* **ITPD bilateral trade data**: "data\source_data\itpd_e_r01\ITPD_E_R01.csv"
* **ITPD sector concordance (ITPD to FAO FCL)**: "data\source_data\itpd_e_r01\ITPD Agriculture FCL concordance.txt"
* **FAO to HS concordance**: "data\source_data\\FAOSTAT_expanded_item_codes_8-19-2021.csv"

MFN Tariffs
* **MFN Data (zipped folder of data files)**: "data\source_data\TRAINS Bulk MFN 2005-2019"
* **MFN tariff country codes**: "data\source_data\TRAINS Bulk MFN 2005-2019\\00_TRAINS_country_codes.csv"

NTM Information
* **NTM Data**: 'data\source_data/ntm_hs6_2010_2018v.12/NTM_hs6_2010_2018 v.12.csv'

Country Information
* **EUN membership**: "data\source_data\EUN_member_codes.csv"
* **DGD gravity files**: 
    * "data/source_data/dgd_2.1/release_2.1_2005_2009.csv"
    * "data/source_data/dgd_2.1/release_2.1_2010_2014.csv"
    * "data/source_data/dgd_2.1/release_2.1_2015_2019.csv"
* **GDP Data (more recent than gravity data)**: "data\source_data\WDI GDP\WDI GDP data 08-05-2021.csv"

#### Output

* **List of Countries with NTM data**: analysis/constructed_data/ntm_data_countries.csv (Available at https://doi.org/10.7910/DVN/AGG5T9)
* **Panel of Unilateral data for stage two regressions**: analysis/constructed_data/unilateral_data.csv (Available at https://doi.org/10.7910/DVN/AGG5T9)
* **Bilateral gravity panel for first stage**: analysis/constructed_data/itpd_bilateral_panel.csv (Available at https://doi.org/10.7910/DVN/AGG5T9)


---
### Step 2: Estimate first stage gravity Models
* **Script**: analysis/gravity_analysis.py

This script prepares a collection of gravity samples and estimates, sector-by-sector, six different gravity specifications. It also produces the table of the distribution of NTM coverage.

Given the number of models that must be estimated, the process can be significantly sped up by dividing the analysis into different groups of sectors and estimating simultaneously in different Python instances. Visual code, which can open up a bunch of separate terminals, worked well for this.

#### Data Inputs
Data inputs all created from Step 1.
* **Bilateral gravity panel**: analysis/constructed_data/itpd_bilateral_panel.csv (Created in Step 1 or available at https://doi.org/10.7910/DVN/AGG5T9)
* **Unilateral data**:  analysis/constructed_data/unilateral_data.csv (Created in Step 1 or available at https://doi.org/10.7910/DVN/AGG5T9)
* **List of countries with NTM data**: analysis\constructed_data\\ntm_data_countries.csv (Created in Step 1 or available at https://doi.org/10.7910/DVN/AGG5T9)

#### Other Inputs
* **root_directory**: Specify root directory for project, which allows for module imports from other scripts. (```root_directory = D:\work\Peter_Herman\projects\\ntm_ave_analysis\\```)
* **TraderRanking.py**: Import TraderRanking function from economic_analysis_tools\\data_analysis\\TraderRanking.py" 
* **Stata installation**: Specify directory of Stata 17 executable (```stata_setup.config("C:\\Program Files\\Stata17", "mp")```)

#### Output
The analysis completes a series of 6 sets of gravity specifications. Each specification produces two files for each sector (three in the HLY specification).

* **<sector_id>_fe_ests.csv**:  A csv file containing the importer and exporter fixed effect estimates for each country (stored as a bilateral panel)
* **<sector_id>_param_ests.dta**: A table of regression estimates, std. errs., etc. for each sector. Stored as a Stata table so that the estimates for both Stata and Python produced estimates match in format.
* (DT2_HLY only) **<sector_id>_log_with_delta_se.log**: A Stata log file of the regression that includes a printout of the standard errors calculated from via delta method. 

The 6 different gravity specifications store their output in different directories within "analysis/results".
* **DA2_expend**: Domestic Trade, all countries, 2-way gravity, trade/expenditures
* **DT2_expend**: Domestic Trade, top countries, 2-way gravity, trade/expenditures
* **DT2_gdp**: Domestic Trade, top countries, 2-way gravity, trade/gdp
* **DT2_hly**: Domestic Trade, top countries, 2-way gravity, ntm and tariff covariates 
* **DT3_expend**: Domestic Trade, top countries, 3-way gravity, trade/expenditures
* **FT2_expend**: Foreign Trade, top countries, 3-way gravity, trade/expenditures

Finally, it outputs a table summarizing the presence of NTMs across sectors in different years.
* **NTM Coverage summary**: analysis/results/summary_ntm_coverage.csv

---
### Step 3: Estimate second stage models
* **Script:** analysis/ave_analysis.py

This script performs most of the analysis presented in the paper, including two-stage estimates, HLY estimates, total NTM costs, aggregate trade cost estimates, and IV estimates.

#### Data Inputs
* **Gravity model estimates**: (Created in Step 2)
    * analysis/results/DA2_expend
    * analysis/results/DT2_expend
    * analysis/results/DT2_gdp
    * analysis/results/DT2_HLY
    * analysis/results/DT3_expend
    * analysis/results/FT2_expend
* **Bilateral gravity panel**: analysis/constructed_data/itpd_bilateral_panel.csv (Created in Step 1 or available at https://doi.org/10.7910/DVN/AGG5T9)
* **Unilateral data**:  analysis/constructed_data/unilateral_data.csv (Created in Step 1 or available at https://doi.org/10.7910/DVN/AGG5T9)
* **List of countries with NTM data**: analysis\constructed_data\\ntm_data_countries.csv (Created in step 1 or available at https://doi.org/10.7910/DVN/AGG5T9)
* **Source GDP Data**: "data\source_data\WDI GDP\WDI GDP data 08-05-2021.csv" 
* **Sector labels**: = "data/source_data/itpd_e_r01/ag_sector_labels.csv"


#### Other Inputs
* **root_directory**: Specify root directory for project, which allows for module imports from other scripts. (```root_directory = D:\work\Peter_Herman\projects\\ntm_ave_analysis\\```)
* **format_regression_table.py**: Import TraderRanking function from economic_analysis_tools\\data_analysis\\format_regression_table.py" 
* **Stata installation**: Specify directory of Stata 17 executable (```stata_setup.config("C:\\Program Files\\Stata17", "mp")```)

#### Output
* **Table of all second stage estimates**: analysis/results/preferred_regressions.<csv/tex> 
* **Table of DT2 estimates for paper**: analysis/results/DT2_table.tex
* **Table of common second stage applied to different gravity specifications for paper**: analysis/results/multi_grav_table.tex
* **Total NTM estimates summarized by sector**: 
    * (Unformatted) analysis/results/total_ntm_summary_by_industry.csv
    * (Formatted) analysis/results/total_ntm_summary_by_industry.tex
* **HLY/Two-stage comparisons**: 
    * (Unformatted) analysis/results/2_stage_HLY_comparison.csv
    * (Formatted) analysis/results/2_stage_HLY_comparison.tex
* **Summary of aggregate cost estimates**: analysis/results/ave_summary_stats.tex
* **IV estimates and regression log**:
    * analysis/results/iv_ests_ntm_mfn.csv
    * analysis/results/iv_ests_ntm_mfn.log


