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
import sys, glob, os, subprocess
import pandas as pd
import numpy as np
import scipy.stats as ss
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineMcf10aCycif as P

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
rSource = 'pipeline/scripts/pipeline-mcf10a-cycif.R'
r.source(rSource)

##### 3. Files #####
concatenatedExpressionFile = 's1-data.dir/HMS_Datasets-merged_filtered.txt'
differentialExpressionFiles = glob.glob('s2-gene_differential_expression.dir/*')

#######################################################
#######################################################
########## S1. Process Data
#######################################################
#######################################################

#############################################
########## 1. Process dataset
#############################################

@follows(mkdir('s1-data.dir/ungrouped'))

@transform(glob.glob('../rawdata/cycif/*/*_data_file.xlsx'),
		   regex(r'../rawdata/cycif/(.*)/.*'),
		   add_inputs(r'../rawdata/cycif/\1/\1_antibody_LUT.xlsx'),
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
	   concatenatedExpressionFile)

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
########## S2. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Gene Differential Expression
#############################################

def deJobs():

	# Get genes
	colnames = subprocess.check_output("head -n1 {concatenatedExpressionFile}".format(**globals()), shell=True).strip().split('\t')
	genes = [x for x in colnames if x not in ['DrugName', 'Conc', 'well', 'fraction', 'timepoint'] ]

	# Generate outfiles
	for gene in genes:
		yield [concatenatedExpressionFile, 's2-gene_differential_expression.dir/{gene}-differential_expression.txt'.format(**locals())]

@follows(mkdir('s2-gene_differential_expression.dir'))

@files(deJobs)

def getGeneDifferentialExpression(infile, outfile):

	# Get gene symbol
	gene_symbol = os.path.basename(outfile)[:-len('-differential_expression.txt')]

	# Get expression data
	expression_dataframe = pd.read_table(infile)[['DrugName', 'Conc', 'fraction', 'timepoint', gene_symbol]]
	expression_dataframe[gene_symbol] = np.log10(expression_dataframe[gene_symbol]+1)

	# Get conditions
	condition_dataframe = expression_dataframe.drop(gene_symbol, axis=1).drop_duplicates().query('DrugName!="DMSO"').reset_index(drop=True)

	# Set index
	expression_dataframe.set_index(['DrugName', 'Conc', 'fraction', 'timepoint'], inplace=True)

	# Get results
	results_list = []

	# Loop through conditions
	for drug_name, concentration, fraction, timepoint in condition_dataframe.as_matrix():

		# Get values to compare
		treated = expression_dataframe.loc[drug_name, concentration, fraction, timepoint][gene_symbol]
		untreated = expression_dataframe.loc['DMSO', 0, fraction, timepoint][gene_symbol]

		# Run test
		test_results = ss.mannwhitneyu(x = treated, y = untreated)

		# LogFC
		logfc = np.log2(np.median([10**x for x in treated])/np.median([10**x for x in untreated]))

		# Add results
		results_list.append({'drug_name': drug_name, 'concentration': concentration, 'fraction': fraction, 'timepoint': timepoint, 'pvalue': test_results.pvalue, 'statistic': test_results.statistic, 'logfc': logfc})

	# Convert to dataframe
	result_dataframe = pd.DataFrame(results_list)

	# Fix
	result_dataframe['bonferroni_pvalue'] = [x*len(result_dataframe.index) if x*len(result_dataframe.index) < 1 else 1 for x in result_dataframe['pvalue']]

	# Write
	result_dataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S3. Signature Matrix
#######################################################
#######################################################

#############################################
########## 1. Signature Matrix
#############################################

@follows(mkdir('s3-signature_matrix.dir'))

# def signatureMatrixJobs():
# 	infiles = differentialExpressionFiles
# 	for fraction in ['cytosolic', 'nuclear']:
# 		for timepoint in ['24h', '48h', '72h']:
# 			outfile = 's3-signature_matrix.dir/cycif-{fraction}_{timepoint}-signature_matrix.txt'.format(**locals())
# 			yield [infiles, outfile]

# @files(signatureMatrixJobs)

def buildSignatureMatrix(infiles, outfile):

# Read data
merged_dataframe = pd.concat([pd.read_table(x) for x in infiles])

# Fix drug labels
merged_dataframe['drug_name'] = [x.lower() for x in merged_dataframe['drug_name']]

# Get conditions
fraction, timepoint = os.path.basename(outfile).split('-')[1].split('_')

# Filter
filtered_dataframe = merged_dataframe.query('timepoint == "{timepoint}" and fraction == "{fraction}"'.format(**locals()))

# Get signature column
filtered_dataframe['signature'] = ['_'.join([x, y, z]) for x, y, z in filtered_dataframe[['signature', 'signature', 'signature']].as_matrix()]

# Filter

	print infiles, outfile


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
