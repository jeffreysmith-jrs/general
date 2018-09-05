#This function calculates abundance based chao similarity as a metric for numpy pdist function
#Created by Jeffrey R. Smith
#Last edited September 5, 2018


#Import necessary library
import numpy as np


def chao(j, k):
	
	#Find species shared between the two sites
	shared	=	np.argwhere((j>=1) & (k>=1)).flatten()
	
	#If there are no shared species return dissimilarity of 1
	if len(shared > 0):
		
		
		#Calculate the total number of individuals in site J and K of all species
		Nj	=	np.sum(j)
		Nk	=	np.sum(k)
		
		#Calculate the number of individuals in site J and K of species in both sites
		Cj	=	np.sum(j[shared])
		Ck	=	np.sum(k[shared])
		

		#Find speices that only occur with 1 individual in site in either J or K that are also present in the other site		
		jones	=	np.argwhere((j>=1) & (k==1)).flatten()
		kones	=	np.argwhere((j==1) & (k>=1)).flatten()
		
		#Find speices that only occur with 2 individuals in site in either J or K that are also present in the other site		
		jtwos	=	np.argwhere((j>=1) & (k==2)).flatten()
		ktwos	=	np.argwhere((j==2) & (k>=1)).flatten()
		
		#Sum up the number of singletons
		a1j	=	len(jones)
		a1k	=	len(kones)
		
		#Sum up doubletons (set to one if there are none as specified in the original paper)
		a2j	=	max(len(jtwos), 1)
		a2k	=	max(len(ktwos), 1)
		
		#Calculate the total number of individuals in site J of species that occur with 1 individual in site K and vice versa
		Sj	=	np.sum(j[jones])
		Sk	=	np.sum(k[kones])
		
		#Calculate indices
		U	=	min(Cj/Nj + (Nk -1)/Nk * a1j/(2*a2j) * Sj/Nj, 1)
		V	=	min(Ck/Nk + (Nj -1)/Nj * a1k/(2*a2k) * Sk/Nk, 1)	

		z = U*V / (U + V - U*V)
		return 1-z

	else:
		return 1
		
#Sample command to run this script
#statDistances	=	scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(communityMatrix, chao))
