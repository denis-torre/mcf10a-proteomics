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
from ruffus import *
import sys, glob, os
import pandas as pd
import numpy as np
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineMcf10aProteomics as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
dataset_names = {
	'HMS_Dataset_20303': 'cytosolic-24h',
	'HMS_Dataset_20304': 'nuclear-24h',
	'HMS_Dataset_20305': 'cytosolic-48h',
	'HMS_Dataset_20306': 'nuclear-48h',
	'HMS_Dataset_20307': 'cytosolic-72h',
	'HMS_Dataset_20308': 'nuclear-72h',
}

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-mcf10a-proteomics.R'
r.source(rSource)

##### 3. Files #####

#######################################################
#######################################################
########## S1. Process Data
#######################################################
#######################################################

#############################################
########## 1. Process dataset
#############################################

@follows(mkdir('s1-data.dir/ungrouped'))

@transform(glob.glob('rawdata.dir/*/*_data_file.xlsx'),
		   regex(r'rawdata.dir/(.*)/.*'),
		   add_inputs(r'rawdata.dir/\1/\1_antibody_LUT.xlsx'),
	       r's1-data.dir/ungrouped/\1_filtered.txt')

def processDatasets(infiles, outfile):

	# Split infiles
	expression_infile, antibody_infile = infiles

	# Get antibody dataframe
	antibody_dataframe = pd.read_excel(antibody_infile, index_col='Data Column')

	# Define antibody dict
	antibody_dict = {}

	# Loop through rows
	for dataColumn in antibody_dataframe.index:

		# If has name
		if type(antibody_dataframe.loc[dataColumn, 'Protein Readout Name(s)']) == unicode:
			antibody_dict[dataColumn] = antibody_dataframe.loc[dataColumn, 'Protein Readout Name(s)'].split(' ')[0]
		elif 'DAPI' in dataColumn:
			antibody_dict[dataColumn] = dataColumn
		else:
			antibody_dict[dataColumn] = antibody_dataframe.loc[dataColumn, 'Other Reagent Name']

	# Read data
	cycif_dataframe = pd.read_excel(expression_infile).set_index(['DrugName', 'Conc', 'well'])[antibody_dict.keys()].rename(columns=antibody_dict)

	# Save data
	cycif_dataframe.to_csv(outfile, sep='\t')	

#############################################
########## 2. Concatenate datasets
#############################################

@merge(processDatasets,
	   's1-data.dir/HMS_Datasets-merged_filtered.txt')

def concatenateDatasets(infiles, outfile):

	# Initialize list
	dataframe_list = []

	# Loop through infiles
	for infile in infiles:

		# Read data
		cycif_dataframe = pd.read_table(infile)

		# Get dataset name
		dataset_name = os.path.basename(infile)[:-len('_filtered.txt')]

		# Get dataset info
		fraction, timepoint = dataset_names[dataset_name].split('-')

		# Add info
		cycif_dataframe['fraction'] = fraction
		cycif_dataframe['timepoint'] = timepoint

		# Append
		dataframe_list.append(cycif_dataframe)

	# Concantate
	concatenated_dataframe = pd.concat(dataframe_list)

	# Save
	concatenated_dataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S1. Get coexpression
#######################################################
#######################################################

#############################################
########## 1. Get coexpression
#############################################

@follows(mkdir('s2-coexpression.dir/ungrouped'))

@transform(processDatasets,
		   regex(r'.*/(.*)_filtered.txt'),
		   r's2-coexpression.dir/ungrouped/\1_coexpression.txt')

def getGeneCoexpression(infile, outfile):

	# Read data
	cycif_dataframe = pd.read_table(infile, index_col='DrugName')

	# Z-score transform
	zscore_dataframe = (cycif_dataframe - cycif_dataframe.mean())/cycif_dataframe.std()

	# Create empty list
	dataframe_list = []

	# Loop through drugs
	for drugName in zscore_dataframe.index.unique():
	    
	    # Get subset
	    drug_dataframe = zscore_dataframe.loc[drugName].T
	    
	    # Get pairwise correlations
	    correlation_dataframe = P.getPairwiseGeneCorrelations(drug_dataframe)
	    
	    # Add drug name
	    correlation_dataframe['drug_name'] = drugName
	    
	    # Append
	    dataframe_list.append(correlation_dataframe)
	    
	# Concatenate
	all_correlation_dataframe = pd.concat(dataframe_list)

	# Save
	all_correlation_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Merge coexpression
#############################################

@merge(getGeneCoexpression,
	   's2-coexpression.dir/coexpression-concatenated.txt')

def mergeGeneCoexpression(infiles, outfile):
 
	# Define dataframe list
	dataframe_list = []

	# Loop through infiles
	for infile in infiles:

		# Read dataframe
		correlation_dataframe = pd.read_table(infile, index_col='drug_name')

		# Subtract DMSO
		# subtracted_dataframe = P.subtractDMSO(correlation_dataframe)

		# Get dataset name
		dataset_name = os.path.basename(infile)[:-len('_coexpression.txt')]

		# Get dataset info
		fraction, timepoint = dataset_names[dataset_name].split('-')

		# Add info
		correlation_dataframe['fraction'] = fraction
		correlation_dataframe['timepoint'] = timepoint

		# Append
		dataframe_list.append(correlation_dataframe)

	# Concatenate
	concatenated_dataframe = pd.concat(dataframe_list)

	# Save
	concatenated_dataframe.to_csv(outfile, sep='\t', index=True)

#############################################
########## 3. Cast by drug
#############################################

@transform(mergeGeneCoexpression,
		   suffix('concatenated.txt'),
		   'by_drug.txt')

def castByDrug(infile, outfile):

	# Read infile
	concatenated_dataframe = pd.read_table(infile)

	# Add new column
	concatenated_dataframe['condition'] = concatenated_dataframe['source']+'__'+concatenated_dataframe['target']+'__'+concatenated_dataframe['drug_name']+'__'+concatenated_dataframe['fraction']

	# Cast
	cast_dataframe = concatenated_dataframe.pivot(index='condition', columns='timepoint', values='spearman_r')

	# Sort values
	sorted_rows = cast_dataframe.apply(np.var, 1).sort_values(ascending=False).index

	# Sort dataframe
	cast_dataframe = cast_dataframe.loc[sorted_rows]

	# # Split conditions
	cast_dataframe.insert(0, 'fraction', [x.split('__')[3] for x in cast_dataframe.index])
	cast_dataframe.insert(0, 'drug_name', [x.split('__')[2] for x in cast_dataframe.index])
	cast_dataframe.insert(0, 'target', [x.split('__')[1] for x in cast_dataframe.index])
	cast_dataframe.insert(0, 'source', [x.split('__')[0] for x in cast_dataframe.index])

	# Write
	cast_dataframe.to_csv(outfile, sep='\t', index=False)


#############################################
# ########## 4. Get difference
# #############################################

# @transform(mergeGeneCoexpression,
# 		   suffix('concatenated.txt'),
# 		   'by_drug.txt')

# def castByDrug(infile, outfile):

# 	# Read infile
# 	concatenated_dataframe = pd.read_table(infile)

# 	# Add new column
# 	concatenated_dataframe['condition'] = concatenated_dataframe['source']+'__'+concatenated_dataframe['target']+'__'+concatenated_dataframe['drug_name']+'__'+concatenated_dataframe['fraction']

# 	# Cast
# 	cast_dataframe = concatenated_dataframe.pivot(index='condition', columns='timepoint', values='pearson_r')

# 	# Sort values
# 	sorted_rows = cast_dataframe.apply(np.var, 1).sort_values(ascending=False).index

# 	# Sort dataframe
# 	cast_dataframe = cast_dataframe.loc[sorted_rows]

# 	# # Split conditions
# 	cast_dataframe.insert(0, 'fraction', [x.split('__')[3] for x in cast_dataframe.index])
# 	cast_dataframe.insert(0, 'drug_name', [x.split('__')[2] for x in cast_dataframe.index])
# 	cast_dataframe.insert(0, 'target', [x.split('__')[1] for x in cast_dataframe.index])
# 	cast_dataframe.insert(0, 'source', [x.split('__')[0] for x in cast_dataframe.index])

# 	# Write
# 	cast_dataframe.to_csv(outfile, sep='\t', index=False)


#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=4, verbose=1)
print('Done!')
