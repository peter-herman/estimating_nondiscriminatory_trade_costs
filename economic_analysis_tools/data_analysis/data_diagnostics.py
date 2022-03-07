__Author__ = "Peter Herman"
__Project__ = "ntm_ave_analysis"
__Created__ = "November 22, 2021"
__Description__ = ''' '''



class CompareIdentifiers(object):
    def __init__(self,
                 dataframe_a = None,
                 dataframe_b = None,
                 code_columns_a = [],
                 code_columns_b = []):
        '''
        A class that compares the identifiers in two datasets to help diagnose merge quality
        Args:
            dataframe_a: (Pandas DataFrame)
                A dataframe containing one set of identifiers
            dataframe_b: (Pandas DataFrame)
                A dataframe containing another set of identifiers
            code_columns_a: (List[str])
                A list of a column name or names containing identifiers (identifiers in multiple columns are combined)
            code_columns_b: (List[str])
                A list of a column name or names containing identifiers (identifiers in multiple columns are combined)
        Attributes:
            codes_a: (List[str])
                A list of identifier codes in dataframe_a
            codes_b: (List[str])
                A list of identifier codes in dataframe_b
            a_not_in_: (List[str])
                A list of identifier codes in dataframe_a but not in dataframe_b
            b_not_in_: (List[str])
                A list of identifier codes in dataframe_b but not in dataframe_a
            in_both: (List[str])
                A list of identifier codes in both dataframe_a and dataframe_b (i.e. the intersection)
            in_either: (List[str])
                A list of identifier codes in either dataframe_a and dataframe_b (i.e. the union)
            code_merge: (Pandas DataFrame)
                A DataFrame with the codes from both matched together
            unmatched_a: (Pandas DataFrame)
                Codes in a that are not matched with b (i.e. codes_b is nan)
            unmatched_b: (Pandas DataFrame)
                Codes in b that are not matched with a (i.e. codes_a is nan)
            unmatched: (Pandas DataFrame)
                All codes in both dataframes that are unmatched
        Methods:
            summary(self):
                Prints basic summary information
        '''
        if not isinstance(code_columns_a, list):
            code_columns_a = [code_columns_a]
        if not isinstance(code_columns_b, list):
            code_columns_b = [code_columns_b]

        codes_a = set()
        for col in code_columns_a:
            temp = dataframe_a[col].unique().tolist()
            codes_a = codes_a.union(set(temp))

        self.codes_a = list(codes_a)

        codes_b = set()
        for col in code_columns_b:
            temp = dataframe_b[col].unique().tolist()
            codes_b = codes_b.union(set(temp))

        self.codes_b = list(codes_b)

        self.a_not_in_b = list(codes_a - codes_b)
        self.b_not_in_a = list(codes_b - codes_a)
        self.in_both = list(codes_a.intersection(codes_b))
        self.in_either = list(codes_a.union(codes_b))

        merge_a = pd.DataFrame(list(codes_a), columns = ['code_a'])

        merge_b = pd.DataFrame(list(codes_b), columns = ['code_b'])


        self.code_merge = merge_a.merge(merge_b, how = 'outer', left_on = ['code_a'], right_on=['code_b'])
        self.unmatched_a = self.code_merge.loc[self.code_merge['code_b'].isnull(),:]
        self.unmatched_b = self.code_merge.loc[self.code_merge['code_a'].isnull(), :]
        self.unmatched = pd.concat([self.unmatched_a, self.unmatched_b])

        self._summary_text = [("Number of 'a' codes: " + str(len(self.codes_a))),
                              ("Number of 'b' codes: " + str(len(self.codes_b))),
                              ("Number of codes in both: " + str(len(self.in_both))),
                              ("Number of codes in 'a' but not 'b': " + str(len(self.a_not_in_b))),
                              ("Number of codes in 'b' but not 'a': " + str(len(self.b_not_in_a)))]
        for text in self._summary_text:
            print(text)

    def summary(self):
        for text in self._summary_text:
            print(text)
        return None

    def __repr__(self):
        strg = self._summary_text
        return "{} \n" \
               "{} \n" \
               "{} \n" \
               "{} \n" \
               "{} ".format(strg[0],strg[1],strg[2],strg[3],strg[4])