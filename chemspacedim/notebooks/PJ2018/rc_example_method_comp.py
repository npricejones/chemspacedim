from spectralspace.analysis.empca_residuals import *

# First and second run - uncomment below, comment out bmask = 7935
#bmask = 4351
# Third run - uncomment below
bmask = 7935
datadir = './data'
maxsamp = 5

redclump = empca_residuals('apogee','red_clump',maskFilter,ask=True,degree=2,badcombpixmask=bmask,datadir = datadir)

redclump.samplesplit(fullsamp=False,subsamples=25,maxsamp=maxsamp,varfuncs=[np.ma.var,meanMed],ctmnorm=None)
