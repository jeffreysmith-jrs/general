#Reads in a Maxent summary file and generates the SDM graph (working-ish)
#Jeffrey Smith
#August 21,2018

#Import libraries
import os, codecs
import numpy as np
import pandas as pd
from glob import glob
from math import exp


def postMaxent(outputDirectory, species):
	if isinstance(species, basestring):
		species	=	[species]
	
	for t in range(len(species)):
		species	=	str(species[t])


		#Read in maxent HTML
		htmlFile		=	codecs.open(str(outputDirectory + species + '.html'), 'r')
		maxentHTML		=	 htmlFile.read()
		htmlFile.close()

		#Get the variables used in the Maxent Run
		batchLineCommand	=	maxentHTML.split('\n')[-2]
		batchLineCommand	=	batchLineCommand.replace('<br>', ' ') #End on a space


		#Figure out where in the string things start and extract relevant information
		startSamples 		=	batchLineCommand.find('samplesfile') + 12
		endSamples			=	batchLineCommand[startSamples:].find(' ') + startSamples
		samplesFile			=	batchLineCommand[startSamples:endSamples]

		#Get environmental layers used to make the maxent model
		startEnv			=	batchLineCommand.find('environmentallayers') + 20 
		endEnv				=	batchLineCommand[startEnv:].find(' ') + startEnv
		inputAsciis			=	batchLineCommand[startEnv:endEnv]

		#Deal with replicated runs
		startReplicates		=	batchLineCommand.find('replicates') + 11 
		if startReplicates != 10:
		
			#This became an issue if replicates was the last thing listed in the command line
			#So I had to add to search for either the next space or the line break
			endReplicates		=	batchLineCommand[startReplicates:].find(' ')
			if endReplicates > 0:
				endReplicates	=	endReplicates + startReplicates
			else:
				endReplicates		=	batchLineCommand[startReplicates:].find(' ') + startReplicates


			
			replicates			=	int(batchLineCommand[startReplicates:endReplicates])
		
		#Unreplicated runs
		else:
			replicates = 1

		#Figure out which, if any asciis were excluded from analysis
		toExclude			=	batchLineCommand
		excludedLayers		=	[]
		while toExclude.find('-N') > 0:
			startExclude	=	toExclude.find('-N') + 3 
			endExclude		=	toExclude[startExclude:].find(' ') + startExclude
			excludeLayer	=	toExclude[startExclude:endExclude]
			toExclude		=	toExclude[endExclude:]
			excludedLayers	=	excludedLayers + [excludeLayer]

		#Read in dataframes	
		#Make blank dictionary to hold them all
		predictors				=	{}
		environmentalLayers		=	glob(str(inputAsciis + '\\*.asc'))
		for i in range(len(environmentalLayers)):
			asciiName	=	environmentalLayers[i].split('\\')[-1].split('.')[0]
			
			#Make sure those asciis were not excludled from the analysis
			if asciiName not in excludedLayers:
				predictors[asciiName]	=	np.genfromtxt(environmentalLayers[i], skip_header = 6)
				predictors[asciiName][predictors[asciiName] == -9999] = np.nan

				

		replicatedRuns	=	np.zeros(shape = (replicates, len(predictors[asciiName]), len(predictors[asciiName][0])))
		for k in range(replicates):
					
			#Make blank basemap
			fx	=	np.zeros_like(predictors[asciiName]) * predictors[asciiName]
			fx	=	fx.astype(np.float64)
					
			#Read in lambdas file from maxent
			
			if replicates > 1:
				fname		=	str(outputDirectory + species + '_' + str(k) + '.lambdas')
			else:
				fname		=	str(outputDirectory + species +  '.lambdas')
			lambdas		=	pd.read_csv(fname, names = ['variable', 'lambda', 'min', 'max'], dtype = {'varaible':str})		
			
			#Extract constants
			linearPredictorNormalizer	=	float(lambdas.loc[lambdas['variable'] == 'linearPredictorNormalizer']['lambda'])
			densityNormalizer			=	float(lambdas.loc[lambdas['variable'] == 'densityNormalizer']['lambda'])
			numBackgroundPoints			=	float(lambdas.loc[lambdas['variable'] == 'numBackgroundPoints']['lambda'])
			entropy						=	float(lambdas.loc[lambdas['variable'] == 'entropy']['lambda'])

			#Drop rows with constants
			lambdas = 	lambdas[:-4]
		
		
		
			#Check for catagorical features:
			catagoricalFeatures	=	lambdas[lambdas['variable'].str.contains('\=')]
			lambdas.drop(catagoricalFeatures.index, inplace = True)
			catagoricalFeatures	=	catagoricalFeatures.reset_index(drop = True)
			
			

			if len(catagoricalFeatures) >= 1:
				for i in range(len(catagoricalFeatures)):
					layer, cutoff		=	catagoricalFeatures['variable'][i].strip('/(').strip('/)').split('=')
					x					=	predictors[layer]
					
					category	=	np.where(x == float(cutoff), catagoricalFeatures['lambda'][i], 0)
					fx 	+= 	category
				
				
				
			#Check for threshold features:
			thresholdFeatures	=	lambdas[lambdas['variable'].str.contains('\<')]
			lambdas.drop(thresholdFeatures.index, inplace = True)
			thresholdFeatures	=	thresholdFeatures.reset_index(drop = True)
			
			

			if len(thresholdFeatures) >= 1:
				for i in range(len(thresholdFeatures)):
					cutoff, layer		=	thresholdFeatures['variable'][i].strip('/(').strip('/)').split('<')
					x					=	predictors[layer]
					
					
					threshold	=	np.where(x < float(cutoff), 0, thresholdFeatures['lambda'][i])
					fx 	+= 	threshold
			
			
			
			
			
			#Check for reverse hinge features
			reverseFeatures	=	lambdas[lambdas['variable'].str.contains('`')]
			lambdas.drop(reverseFeatures.index, inplace = True)
			reverseFeatures	=	reverseFeatures.reset_index(drop = True)
			
			if len(reverseFeatures) >= 1:
				for i in range(len(reverseFeatures)):
					a		=	reverseFeatures['variable'][i].split('`')[1]
					x		=	predictors[a]
					
					hinge	=	np.where(x < reverseFeatures['max'][i], (reverseFeatures['lambda'][i] *
																		(reverseFeatures['max'][i] - x) /
																		(reverseFeatures['max'][i] - reverseFeatures['min'][i])),
																		0)
					fx += hinge	
					
		
			#Check for y 
			forwardFeatures	=	lambdas[lambdas['variable'].str.contains("'")]
			lambdas.drop(forwardFeatures.index, inplace = True)
			forwardFeatures	=	forwardFeatures.reset_index(drop = True)
			
			if len(forwardFeatures) >= 1:
				for i in range(len(forwardFeatures)):
					a		=	forwardFeatures['variable'][i].split("'")[1]
					x		=	predictors[a]
					
					
					hinge	=	np.where(x >= forwardFeatures['min'][i],	(forwardFeatures['lambda'][i] *
																		(x - forwardFeatures['min'][i]) /
																		(forwardFeatures['max'][i] - forwardFeatures['min'][i])),
																		0)
					fx 	+= 	hinge
			

			#Check for multiplicative features
			multiplicativeFeatures	=	lambdas[lambdas['variable'].str.contains('\*')]
			lambdas.drop(multiplicativeFeatures.index, inplace = True)
			multiplicativeFeatures	=	multiplicativeFeatures.reset_index(drop = True)

			#Calculate contribution of multiplicative features
			if len(multiplicativeFeatures) >= 1:
				for i in range(len(multiplicativeFeatures)):
					a, b 	=	multiplicativeFeatures['variable'][i].split('*')
					x1		=	predictors[a]
					x2		=	predictors[b]

					fx		+=	(	multiplicativeFeatures['lambda'][i] *
									(x1 * x2 - multiplicativeFeatures['min'][i]) /
									(multiplicativeFeatures['max'][i] - multiplicativeFeatures['min'][i]))

			#Check for squared features
			squareFeatures	=	lambdas[lambdas['variable'].str.contains('\^')]
			lambdas.drop(squareFeatures.index, inplace = True)
			squareFeatures	=	squareFeatures.reset_index(drop = True)

			#Calculate contribution of squared features
			if len(squareFeatures) >= 1:
				for i in range(len(squareFeatures)):
					a		=	squareFeatures['variable'][i].split('^')[0]
					x		=	predictors[a]

					fx		+=	(	squareFeatures['lambda'][i] *	
									(x * x - squareFeatures['min'][i]) / 
									(squareFeatures['max'][i] - squareFeatures['min'][i]))
					
			#Finally do linear features
			lambdas	=	lambdas.reset_index(drop = True)
			if len(lambdas) >= 1:
				for i in range(len(lambdas)):
					a		=	lambdas['variable'][i]
					x		=	predictors[a]

					fx		+=	(	lambdas['lambda'][i] *	
									(x - lambdas['min'][i]) / 
									(lambdas['max'][i] - lambdas['min'][i]))
					
			s	=	fx - linearPredictorNormalizer
			qx	=	np.exp(s) / densityNormalizer
			ls	=	(qx * np.exp(entropy) / (1 + qx * np.exp(entropy)))
			
			replicatedRuns[k,:,:]	=	np.copy(ls)
			
	return replicatedRuns

