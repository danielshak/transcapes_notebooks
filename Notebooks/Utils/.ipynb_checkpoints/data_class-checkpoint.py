import pandas as pd
import numpy as np
#Data class for fast testing of Gpy

class GPY_data():
    def __init__(self,pdDF,features,cell_lines):
        '''
        Class that takes in dataframe of all the cell conditions and filters for various parameters. Creates sub dataframes by filtering for
        missing gene values, missing protein values and individual chromosome level. Dataframe exists as a 3 column df, to get sub columns of
        the index use:
            #test.rna['AT2_04M_F0'].index.get_level_values(0) #selects a single index for multiindex
        
        self.cell_lines has the names of the different conditions, use that as keys to get data from:
            self.rna - (dictionary) rows of genes of nonzero mRNA levels
            self.rna_protein - (dictinary) rows of genes of nonzero mRNA and protein levels
            self.chrm - (dictionary of dictionary) top level dictionary keys are cell lines, second level down are chrm{i}, 1<=i<=21
        
        Various functions for testing gpy models, and uniform plotting:
            -grid to project back onto chrm position
            -grid 
            
        '''
        assert isinstance(pdDF,pd.DataFrame)
        assert isinstance(features, list)
        assert isinstance(cell_lines, list)
        
        self.data = pdDF.copy() #pandas dataframe
        self.data['Gene Length'] = np.log2(self.data[['Gene Length']]) #log2 of gene length
        self.features = features
        self.cell_lines = cell_lines #single list
        
        self.rna = {} #Dictionary containing dataframe for cell conditions with rows of genes of nonzero mRNA levels
        self.rna_chrm = {} #Dictinoary of data frames filtered for each chromosome, contains all genes with nonzero mRNA
        self.rna_protein = {} #Dictionary containing dataframe for cell conditions with rows of genes of nonzero mRNA and protein levels
        self.rna_protein_chrm = {} #Data frames filtered for each chromosome, contains only nonzero mRNA and protein
        num_chrms = 21 #accounting for both x/y chrms
        
        for cell_condition in cell_lines:
            rna_prot = [cell_condition,cell_condition+'_P'] #mRNA and protein names
            self.rna[cell_condition] = self.data[features+rna_prot][self.data[cell_condition]!=0] #removes genes where mRNA = 0
            self.rna_protein[cell_condition] = self.rna[cell_condition][self.rna[cell_condition][cell_condition+'_P']!=0] #removes genes where mRNA and protein levels are 0
            sub_rna_chrm = {}
            sub_chrm = {}
            for i in range(num_chrms):
                sub_rna_chrm[f'chrm{i+1}'] = self.rna[cell_condition][(self.rna[cell_condition]['AvgChrs']>=i) & (self.rna[cell_condition]['AvgChrs']<i+1)] #filters per chromosome level
                sub_chrm[f'chrm{i+1}'] = self.rna_protein[cell_condition][(self.rna_protein[cell_condition]['AvgChrs']>=i) & (self.rna_protein[cell_condition]['AvgChrs']<i+1)] #filters per chromosome level
            self.rna_chrm[cell_condition] = sub_rna_chrm
            self.rna_protein_chrm[cell_condition] = sub_chrm
                