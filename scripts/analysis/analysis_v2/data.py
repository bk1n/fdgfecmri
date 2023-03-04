import pandas as pd
import numpy as np
from textwrap import dedent
from sklearn import impute

# TODO - save these dataframes (incl. df_filt)!
# TODO - split byCANC
# TODO - start on building models


class Data():
    def __init__(self) -> None:
        self.df = pd.read_csv('./data/quant_allPooled.csv')
        #rows to drop:
            # FEC_SUV
            # FEC_SA, FEC_LA
            # FEC_NTR
            # FEC_STAR
            # FEC_SNSA
            # FDG_SUV_PT
            # ADC_PT
            # ADC_NTR
            
    def stats(self):
        '''Print basic summary statistics for the data (e.g. no. of NaNs etc.)'''
        
        print(dedent(
            f'''
            ---- STATS ----
            
            -- Number of NAs in each column:'''))
        print(self.df.isna().sum(axis=0).sort_values(ascending=False))

        print(f'\n -- Number of non-NAs in each column')
        print(self.df.notna().sum(axis=0).sort_values(ascending=True))
        
        print('\n-- Percentage of NAs')
        print((self.df.isna().sum(axis=0) / len(self.df)).sort_values(ascending=False))

        print(f'\n-- Number of rows with 0 NAs: \n{sum(self.df.isna().sum(axis=1) == 0)}')
        print(f'\n-- Mean number of NAs per row: \n{self.df.isna().sum(axis=1).mean():.2f}')
              
        print('\n-- Column stats')
        print(self.df.describe())
        
        print(f'\n-- Cancer type count: \n{self.df["CANC"].value_counts()}')
        print(f'\n-- Hist count: \n{self.df["HIST"].value_counts()}')
            
    def combined(self):
        self.df = self.df.drop(['FEC_SUV', 'FEC_SA', 'FEC_LA', 'FEC_NTR', 'FEC_STAR', 'FEC_SNSA', 'FEC_SUV_PT', 'FDG_SUV_PT', 'ADC_PT', 'ADC_NTR'], axis=1)
        
    def get_endo(self):
        self.df = self.df[self.df['CANC'] == 'endo']
        
    def get_cer(self):
        self.df = self.df[self.df['CANC'] == 'cer']
        
    def remove_columns_less50_null(self):
        columns = (self.df.isnull().sum(axis=0) / len(self.df))
        columns = list(columns[columns < .5].index)
        self.df = self.df.loc[:,columns]
        
    def remove_null_rows(self):
        self.df = self.df.loc[self.df.isnull().sum(axis=1) == 0,:]
        
    def remove_null_rows_without_fdg_adc(self):
        self.df = self.df.loc[self.df['ADC'].notna() | self.df['FDG_SUV'].notna(),:]
        
    def remove_calculated_variables(self):
        self.df = self.df.drop(['FDG_NTR', 'FDG_STAR', 'FDG_SNSA', 'ADC_NTR'], axis=1)
        
    def knn_impute(self, k):
        knn = impute.KNNImputer(n_neighbors=k)
        impute_df = self.df.drop(['CANC', 'HIST'],axis=1)
        dropped_cols = self.df[['CANC', 'HIST']].reset_index(drop=True)
        colnames = ['CANC', 'HIST'] + list(impute_df.columns)
        impute_df = knn.fit_transform(impute_df)
        impute_df = pd.DataFrame(impute_df)
        
        self.df = pd.concat([dropped_cols, impute_df], axis=1)
        self.df = self.df.set_axis(colnames, axis=1, inplace=False)
              
class FDG(Data):
    def __init__(self) -> None:
        super().__init__()
        self.df = self.df.loc[:, self.df.columns.str.contains('ADC|FDG|HIST|CANC')].drop(['ADC_PT', 'ADC_NTR'], axis=1)
        
class FEC(Data):
    def __init__(self) -> None:
        super().__init__()
        self.df = self.df.loc[:, self.df.columns.str.contains('FEC|HIST|CANC')]
    
class ADC(Data):
    def __init__(self) -> None:
        super().__init__()
        self.df = self.df.loc[:, self.df.columns.str.contains('ADC|HIST|CANC')]
        
if __name__ == '__main__':
    # combined = Data()
    # combined.combined()
    # combined.stats()
    
    # fdg = FDG()
    # fdg.get_endo()
    # fdg.stats()
    
    # fdg.remove_null_rows_without_fdg_adc()
    # fdg.knn_impute(5)
    # fdg.stats()
    
    # fec = FEC()
    # fec.get_endo()
    # fec.stats()
    # fec.remove_null_rows()
    # fec.stats()
    
    # adc = ADC()
    # adc.stats()
    # adc.remove_null_rows()
    # adc.stats()
    
    data = Data()
    data.stats()
    # data.remove_null_rows()
    
    data.remove_null_rows_without_fdg_adc()
    data.remove_columns_less50_null()
    data.stats()
    
    # use this dataset for model 
    # first, remove cervical cancer - not enough data
    # then, split the data into train and test splits (re-write code above)
    # impute this data (should give us 191*train_test split training samples total)
