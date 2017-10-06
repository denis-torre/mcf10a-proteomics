#################################################################
#################################################################
###############  ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
import sys
import scipy.stats as ss
import pandas as pd

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
import Support as S 

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

#######################################################
#######################################################
########## S1. Correlation Analysis
#######################################################
#######################################################

#############################################
########## 1. Pairwise Gene Correlations
#############################################

def getPairwiseGeneCorrelations(expression_dataframe):
    
    # Create list
    correlation_list = []
    
    # Get number of genes
    n_genes = len(expression_dataframe.index)

    # Get index of first gene
    for i in range(n_genes):
        
        # Get index of second gene
        for j in range(i+1, n_genes):
            
            # Get x and y
            x = expression_dataframe.iloc[i]
            y = expression_dataframe.iloc[j]
            
            # Get gene names
            source = expression_dataframe.index[i]
            target = expression_dataframe.index[j]
            
            # Index genes
            cor_results = ss.spearmanr(x, y)
            
            # Add
            correlation_list.append({'source': source, 'target':target, 'spearman_r': cor_results[0], 'pvalue': cor_results[1]})
    
    # Convert to dataframe
    correlation_dataframe = pd.DataFrame(correlation_list)
    
    # Return
    return correlation_dataframe

#############################################
########## 2. Subtract DMSO
#############################################

def subtractDMSO(correlation_dataframe):

	# Get DMSO
	dmso_dataframe = correlation_dataframe.loc['DMSO'].set_index(['source', 'target'])

	# Drop DMSO
	correlation_dataframe.drop('DMSO', inplace=True)

	# Initialize list
	dataframe_list = []

	# Loop through drugs
	for drugName in correlation_dataframe.index.unique():
	    
	    # Get subset
	    drug_dataframe = correlation_dataframe.loc[drugName].set_index(['source', 'target'])
	    
	    # Subtract
	    subtracted_dataframe = dmso_dataframe.loc[drug_dataframe.index] - drug_dataframe
	    
	    # Add drug
	    subtracted_dataframe['drug_name'] = drugName
	    
	    # Append
	    dataframe_list.append(subtracted_dataframe)
	    
	# Concatenate
	concatenated_subtracted_dataframe = pd.concat(dataframe_list).reset_index().drop('pvalue', axis=1).rename(columns={'spearman_r': 'spearman_r_diff'})

	# Return
	return concatenated_subtracted_dataframe


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

