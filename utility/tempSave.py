#Function creates a temporary directory and keeps python graphs there
#Use this just like you use plt.savefig() except no need to add a name
#Jeffrey R. Smith
#July 17, 2018

def tempSave(dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None):

	#Import libraries
	import os, time
	from glob import glob
	import matplotlib.pyplot as plt
	import seaborn as sns
	
	#Specify where you want to save your temporary output graphs
	tempDirectory	=	'C:\\Users\\Jeffrey\\Anaconda2\\tempPlots\\'	
	
	#If that directory doesn't exist make it
	if not os.path.exists(tempDirectory):
		os.makedirs(tempDirectory)
	
	#Specify how many files you want to keep at any one time
	numberPlots		=	100	

	#Define extension for temporary save file (default to png if not specified)
	extension	=	str('.' + str(format)) or '.png'
	
	#Save graph with temporary name
	fname		=	str(tempDirectory + str(time.time()) + extension)
	
	#Save figure command
	plt.tight_layout()
	plt.savefig(fname, dpi=dpi, facecolor=facecolor, edgecolor=edgecolor,
        orientation=orientation, papertype=papertype, format=format,
        transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches,
        frameon=frameon)
	
	#If you have more than your specified number of graphs you wanna hang onto, delete the oldest
	currentGraphs			=	glob(str(tempDirectory + '*'))
	if len(currentGraphs)	>	numberPlots:
		os.remove(currentGraphs[0])
		
		
