__Author__ = "Peter Herman"
__Project__ = "misc_tools"
__Created__ = "January 28, 2020"
__Description__ = '''This is a code to pare down a gravity dataset to include only a subset of the most 
    prominent trading countries. This '''

import pandas as pd


class TraderRanking():
    '''
    Determine and rank countries by trade participation by country, year, and sector.

    Args:
        gravity_data: (pd.DataFrame) A gravity dataset containing columns corresponding to importer, exporter, and
            trade values (necessary) and year (optional).
        imp_var_name: (str) The column name containing importer identifiers. Default is 'importer'.
        exp_var_name: (str) The column name containing exporter identifiers. Default is 'exporter'
        trade_var_name: (str) The column name containing trade flows. Default is 'trade_value'.
        year_var_name: (str) The column name containing year identifiers. Default is 'year'.
        sector_var_name: (str) (optional) The column name containing sector identifiers.

    Methods:
        ranking(flow: str = 'both', by_year: bool = False, by_sector: bool = False)
        top_countries(self, flow='both', country_number=None, by_year=False, by_sector=False)
        cumulative_coverage(self, flow='both', cumulative_percentage=None, by_year=False, by_sector = False)

    Returns: A TradeRanking object.

    Examples:
     >>> two_sector_rank = TraderRanking(two_sector_panel, sector_var_name='Item')
     >>> two_sector_rank.ranking(flow = 'imports', by_sector=True, by_year = True)
                  total_imports     share  cumulative_share
        importer
        CHN        2.840445e+08  0.337741          0.337741
        FRA        9.032501e+07  0.107400          0.445142
        USA        5.561909e+07  0.066134          0.511275
        ESP        3.690341e+07  0.043880          0.555155
        TUR        2.788552e+07  0.033157          0.588312
                         ...       ...               ...
        NAM       -1.100453e+05 -0.000131          1.002251
        SYR       -2.345700e+05 -0.000279          1.001972
        AFG       -2.452990e+05 -0.000292          1.001680
        UZB       -6.036890e+05 -0.000718          1.000963
        PER       -8.096125e+05 -0.000963          1.000000
     >>> two_sector_rank.top_countries(flow = 'exports', country_number=8, by_sector=True, by_year=False)
         {'Grapes': ['CHN', 'FRA', 'USA', 'ESP', 'ITA', 'CHL', 'TUR', 'ZAF'],
          'Apples': ['CHN', 'USA', 'IND', 'JPN', 'IRN', 'TUR', 'ITA', 'FRA']}
     >>> two_sector_rank.cumulative_coverage(flow = 'both', cumulative_percentage=0.6, by_sector = True, by_year = True)
        {'Grapes': {'2010': ['FRA', 'CHN', 'ESP', 'USA', 'TUR', 'ITA', 'ZAF'],
          '2011': ['CHN', 'FRA', 'ESP', 'USA', 'TUR'],
          '2012': ['CHN', 'FRA', 'USA', 'ESP', 'ITA', 'TUR', 'ZAF'],
          '2013': ['FRA', 'CHN', 'USA', 'ESP', 'ITA', 'ZAF', 'TUR'],
          '2014': ['CHN', 'FRA', 'ESP', 'USA', 'ITA'],
          '2015': ['CHN', 'FRA', 'USA', 'ESP', 'ITA'],
          '2016': ['FRA', 'CHN', 'USA', 'ESP', 'CHL']},
         'Apples': {'2010': ['CHN'],
          '2011': ['CHN', 'IND', 'USA', 'IRN', 'TUR'],
          '2012': ['CHN', 'USA'],
          '2013': ['CHN', 'USA', 'IRN', 'IND', 'JPN'],
          '2014': ['CHN', 'USA', 'IND', 'JPN'],
          '2015': ['CHN', 'USA', 'IND'],
          '2016': ['CHN', 'IND', 'USA', 'JPN']}}
    '''

    def __init__(self,
                 gravity_data,
                 imp_var_name: str = 'importer',
                 exp_var_name: str = 'exporter',
                 trade_var_name: str = 'trade_value',
                 year_var_name: str = 'year',
                 sector_var_name: str = None):
        self.gravity_data = gravity_data
        self.imp_var_name = imp_var_name
        self.exp_var_name = exp_var_name
        self.trade_var_name = trade_var_name
        self.year_var_name = year_var_name
        self.year_list = self.gravity_data[self.year_var_name].unique().tolist()
        if sector_var_name:
            self.sector_var_name = sector_var_name
            self.sector_list = self.gravity_data[self.sector_var_name].unique().tolist()


    def _get_ranking(self, gravity_data, flow):
        '''
        A private function for constructing a ranking given a supplied dataset or subset and a specified flow type.
        '''
        total_imports = gravity_data.groupby([self.imp_var_name]).agg({self.trade_var_name: 'sum'})
        total_imports.rename(columns={self.trade_var_name: 'total_imports'}, inplace=True)
        total_exports = gravity_data.groupby([self.exp_var_name]).agg({self.trade_var_name: 'sum'})
        total_exports.rename(columns={self.trade_var_name: 'total_exports'}, inplace=True)
        if flow == 'both':
            flow_column = 'total_trade'
            total_trade = pd.merge(left=total_exports, left_index=True,
                                   right=total_imports, right_index=True, how='outer')
            total_trade = total_trade.fillna(0)
            total_trade[flow_column] = total_trade['total_exports'] + total_trade['total_imports']
        if flow == 'exports':
            flow_column = 'total_exports'
            total_trade = total_exports
        if flow == 'imports':
            flow_column = 'total_imports'
            total_trade = total_imports

        total_trade.sort_values([flow_column], ascending=False, inplace=True)
        total_value = sum(total_trade[flow_column])
        total_trade['share'] = total_trade[flow_column] / total_value
        total_trade['cumulative_share'] = total_trade[flow_column].cumsum() / total_value
        total_trade['rank'] = total_trade.reset_index().index + 1
        return total_trade

    def ranking(self, flow: str = 'both', by_year: bool = False, by_sector: bool = False):
        '''
        Generate a DataFrame of countries ranked by the desired trade flow (imports, exports, or both) from largest to
        smallest. The DataFrame also includes data on each countries share of total flows and cumulative share.

        Args:
            flow: (str) Type of flow to base the ranking on ('imports', 'exports', or 'both'). Default is 'both.
            by_year: (bool) If True, ranking is compiled on a year-by-year basis. If False, ranking reflects total flows
                summed across all years. Default is False.
            by_sector: (bool) If True, ranking is compiled on a sector-by-sector basis. If False, ranking reflects total
             flows summed across all sectors. Default is False.

        Returns:
            (DataFrame or Dict[DataFrame] or Dict[Dict[DataFrame]]) A DataFrame or dictionary of DataFrames containing
            the ranking of countries by the selected flow and some additional related info like each countries share of
            the total and the cumulative share. If by_year or by_sector are selected, the return is a dictionary of
            DataFrames keyed by the sector or year IDs. If both are selected, the return is a dictionary keyed by sector
            labels containing dictionaries keyed by year labels with DataFrame attributes.
        '''
        if not by_sector:
            if not by_year:
                return self._get_ranking(self.gravity_data, flow)
            if by_year:
                yearly_dict = dict()
                for yr in self.year_list:
                    year_data = self.gravity_data.loc[self.gravity_data[self.year_var_name] == yr, :].copy()
                    yearly_dict[yr] = self._get_ranking(year_data, flow)
                return yearly_dict
        if by_sector:
            sector_dict = dict()
            for sec in self.sector_list:
                sector_data = self.gravity_data.loc[self.gravity_data[self.sector_var_name] == sec, :].copy()
                if not by_year:
                    sector_dict[sec] = self._get_ranking(sector_data, flow)
                if by_year:
                    yearly_dict = dict()
                    for yr in self.year_list:
                        year_data = sector_data.loc[sector_data[self.year_var_name] == yr, :].copy()
                        yearly_dict[yr] = self._get_ranking(year_data, flow)
                    sector_dict[sec] = yearly_dict
            return sector_dict

    def _top_countries(self, ranking, country_number):
        '''
        A private function for computing and returning the top countries given a supplied ranking.
        '''
        return ranking.index[0:country_number].tolist()

    def top_countries(self, flow:str ='both', country_number:int =None,
                      by_year:bool=False, by_sector:bool=False):
        '''
        Create a list of the N-many top trading countries.
        Args:
            flow: (str) Type of flow to base the ranking on ('imports', 'exports', or 'both'). Default is 'both.
            country_number: (int) The number of countries to include in the list of top countries.
            by_year: (bool) If True, ranking is compiled on a year-by-year basis. If False, ranking reflects total flows
                summed across all years. Default is False.
            by_sector: (bool) If True, ranking is compiled on a sector-by-sector basis. If False, ranking reflects total
             flows summed across all sectors. Default is False.

        Returns: A list of countries, a dictionary of lists of counties, or a dictionary of dictionaries containing a
        list of countries. See .ranking() for a description of the by_year and/or by_sector dictionary nesting.

        '''
        ranking = self.ranking(flow, by_year, by_sector)
        if not by_sector:
            if not by_year:
                return self._top_countries(ranking,country_number)
            if by_year:
                yearly_countries = dict()
                for yr in self.year_list:
                    yearly_countries[yr] = self._top_countries(ranking[yr], country_number)
                return yearly_countries
        if by_sector:
            sector_dict = dict()
            for sec in self.sector_list:
                sec_ranking = ranking[sec]
                if not by_year:
                    sector_dict[sec] = self._top_countries(sec_ranking,country_number)
                if by_year:
                    year_dict = dict()
                    for yr in self.year_list:
                        year_ranking = sec_ranking[yr]
                        year_dict[yr] = self._top_countries(year_ranking, country_number)
                    sector_dict[sec] = year_dict
            return sector_dict


    def _cumulative_coverage(self, ranking, cumulative_percentage):
        '''
        A private function for computing a list of countries from s aupplied ranking that combine to cover the supplied
        cumulative share of trade flows.
        '''
        applicable_rows = ranking.loc[ranking['cumulative_share'] < cumulative_percentage, :]
        number_of_traders = applicable_rows.shape[0] + 1
        return ranking.index[0:number_of_traders].tolist()

    def cumulative_coverage(self, flow:str ='both', cumulative_percentage:float =None,
                            by_year:bool =False, by_sector:bool = False):
        '''
        Create a list of the top countries that represent a (user supplied) combined percentage of trade. (e.g. the
        top countries that represent 0.9 (90%) of trade).
        Args:
            flow: (str) Type of flow to base the ranking on ('imports', 'exports', or 'both'). Default is 'both.
            cumulative_percentage: (float)
            by_year: (bool) If True, ranking is compiled on a year-by-year basis. If False, ranking reflects total flows
                summed across all years. Default is False.
            by_sector: (bool) If True, ranking is compiled on a sector-by-sector basis. If False, ranking reflects total
             flows summed across all sectors. Default is False.

        Returns: A list of countries, a dictionary of lists of counties, or a dictionary of dictionaries containing a
        list of countries. See .ranking() for a description of the by_year and/or by_sector dictionary nesting.

        '''
        ranking = self.ranking(flow, by_year, by_sector)
        if not by_sector:
            if not by_year:
                return self._cumulative_coverage(ranking, cumulative_percentage)
            if by_year:
                yearly_countries = dict()
                for yr in self.year_list:
                    yearly_countries[yr] = self._cumulative_coverage(ranking[yr], cumulative_percentage)
                return yearly_countries
        if by_sector:
            sector_dict = dict()
            for sec in self.sector_list:
                sec_ranking = ranking[sec]
                if not by_year:
                    sector_dict[sec] = self._cumulative_coverage(sec_ranking, cumulative_percentage)
                if by_year:
                    year_dict = dict()
                    for yr in self.year_list:
                        year_ranking = sec_ranking[yr]
                        year_dict[yr] = self._cumulative_coverage(year_ranking, cumulative_percentage)
                    sector_dict[sec] = year_dict
            return sector_dict