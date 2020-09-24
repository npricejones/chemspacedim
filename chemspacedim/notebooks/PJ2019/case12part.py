"""case12part - run the cluster finding      
Usage:                                                                                                         
    case12part [-h] [-n NUMSTR] [-v VOLUME] [-a INDEX]                    

Options:                                                                                                                 
    -h, --help                              Show this screen   
    -n NUMSTR, --numstr NUMSTR              Total number of stars observed [default: 5e4]                                                
    -v VOLUME, --volume VOLUME              Total annulus volume in cubic kpc [default: 300]
    -a INDEX, --index INDEX                 CMF index [default: -2.1]        

Examples:
    python case12part.py -n 5e4 -v 300 -a -2.1 
"""

import docopt
import numpy as np
from case_template import *

arguments = docopt.docopt(__doc__)

# Read in arguments

# Number of stars observed
nstars = int(float(arguments['--numstr']))
# Power law index of cluster mass function
cmfind = np.fabs(float(arguments['--index']))
# Volume from which cluster members could be sampled
volume = float(arguments['--volume'])

casename = '12-V{0}-a{1}'.format(volume,cmfind)

# run parameters                                                                                                             
sample= 'allStar_chemscrub_teffcut_dr14.npy'#'red_giant_teffcut_dr14.npy' # APOGEE sample to draw from            
abundancefac = 1 # scaling factor for abundance noise                          
specfac = 1e-2 # scaling factor for spectra noise                              
centerfac = 1
suff = 'H' # element denominator                                               
metric = 'precomputed' # metric for distances                                  
fullfitkeys = ['TEFF','LOGG'] # keys for the full fit                          
fullfitatms = []
crossfitkeys = []
crossfitatms = [26] # atomic numbers of cross terms                                                                   
spreadchoice = spreads # choose which abudance spreads to employ

# DBSCAN parameters                                                             

#min_samples = np.array([2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,20,30,40,50])
min_samples = np.array([3])
samples = len(min_samples)
eps = np.array([0.12,0.15,0.19])#0.01,0.02,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.15,0.19,0.24,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
min_samples = np.tile(min_samples,len(eps))
eps = np.repeat(eps,samples)

jobs=8

combelem = ['Mg','Al','Si','S','K','Ca','Ni']

case = caserun()

case.makedata(nstars=nstars,sample=sample,abundancefac=abundancefac,volume=volume,
                 spreadchoice=spreadchoice,specfac=specfac,centerfac=centerfac,
                 centerspr=spreads,genfn=choosestruct,clsind=cmfind,
                 fullfitkeys=fullfitkeys,fullfitatms=fullfitatms,
                 crossfitkeys=crossfitkeys,crossfitatms=crossfitatms,
                 phvary=True,fitspec=True,case=casename,usecenters=True,add=True)

start = time.time()

#case.clustering(case.specinfo.spectra,'spec',eps,min_samples,metric='precomputed',n_jobs=jobs,neighbours = 20,normeps=normeps)
#case.clustering(case.abundances,'abun',eps,min_samples,metric='precomputed',n_jobs=jobs,neighbours = 20,normeps=normeps)
#case.clustering((case.abundances.T[abuninds(elem,combelem)]).T,'reda',eps,min_samples,metric='precomputed',n_jobs=jobs,neighbours = 20,normeps=normeps)
#case.gen_abundances(1,tingspr)
#case.clustering(case.abundances,'tabn',eps,min_samples,metric='precomputed',n_jobs=jobs,neighbours = 20,normeps=normeps)
#case.gen_abundances(1,leungspr)
#case.clustering(case.abundances,'labn',eps,min_samples,metric='precomputed',n_jobs=jobs,neighbours = 20,normeps=normeps)
case.reduction(reduct = PCA, n_components=30)
case.partition(2,case.projectspec,'prin30',eps,min_samples,partitionfns=[case.kmeanspart,case.custompart],partitionkwargs = [{'radiscale':1.0},{'size':int(3e4)}],metric='precomputed',n_jobs=jobs,neighbours = 20,normeps=normeps)
#case.clustering(case.projectspec,'prin30',eps,min_samples,metric='precomputed',n_jobs=jobs,neighbours = 20,normeps=normeps)

end = time.time()
case.finish()
print('Finished desired clustering in {0} seconds'.format(end-start))
