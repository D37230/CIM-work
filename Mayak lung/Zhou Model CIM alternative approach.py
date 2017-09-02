
# coding: utf-8
'''
Evaluate fitted values and compute Q and G from Mmodel coavariate files (rather than the ungroup file)
and compute the CIM mstrix from a dataset of interploated dose realization values

The input files are
1) modcovs.csv  -  model covariates with case counts and person years for ungrouped data 
2) dfreps.csv  -  interpolated dose replicatons for ungrouped data
3) mwds16pualtmods.psv - PSV file we focus on the M " F:M sex ratio model (0)

'''



# load some standard libraries
import os
import sys
import time
import numpy as np   # numerical methods
import h5py  # Hdf5 library
import pandas as pd


# Define mygithub and hdf5path directories and all else should be fine

mygithub = 'C:/Users/Dale/Documents/GitHub/'
hdf5path = 'k:/UralsDosimetry/'

# specify dose repolication and modcov file names
intdrepfn = 'dfreps.csv'
mcfname = 'modcovsfix.csv'
psvname = 'mwds16pufullalt.psv' 

projdir = mygithub + 'CIM-work/Mayak lung/'
modinfodir = hdf5path + 'mayak/mwds16lung/'
# rlzdir = hdf5path + 'mayak/mwds2016/'
os.chdir(projdir)
sys.path.append('../pythonutils')

# import utility functions for this applicaiton

from psvdict2 import modeleval, getmodinfo, calc_fvalues, showfittedcases
from doseintfun import setdoseintinfo, doseinterp
from cimfuncs import mkcim, mkcimdr, getCI, dispcimbnds, mkQmat, modcimprmbnds


# open dose replocatom and modcov files

intdrep = np.genfromtxt(modinfodir + intdrepfn, delimiter=',', names=True)
mwcmc = np.genfromtxt(modinfodir + mcfname, delimiter=',', names=True)
print "Ungrouped record count:", len(mwcmc)

# open psave file and get model information
psavef = open(modinfodir + psvname, 'r')
models = getmodinfo(psavef)
for mno in range(len(models)):
    print "model ",mno,models[mno]['title']
psavef.close()



# ## Set model form and analysis subset


'''
    umod = 0 | 3   # male + F:M sex ratio model
    umod = 1 | 4   # female + M"F sex ratio model
    umod = 2 | 5   # male +_ female dose effect model

'''

umod = 0


tnames = ["Baseline", "Smoking", "External dose", "Pu Surrogate" ,"Pu dose"]
mwcungsub = mwcmc

modidx = umod
tnames = ["Baseline", "Smoking", "External dose", "Pu Surrogate" ,"Pu dose"]
    
print len(mwcungsub), 'records in anlaysis subset'


# ## Get model covariates from modcovs data (there must be a better way)

udbid = mwcmc['udbid'].astype(int)
udbid = np.reshape(udbid,(-1,1))

cases = mwcmc['cases'].astype(int)
cases = np.reshape(cases,(-1,1))

pyr = mwcmc['pyr']
pyr = np.reshape(pyr,(-1,1))

male = mwcmc['male'].astype(int)
male = np.reshape(male,(-1,1))

female = mwcmc['female'].astype(int)
female = np.reshape(female,(-1,1))

con = male + female

lage70 = mwcmc['lage70']
lage70 = np.reshape(lage70,(-1,1))

lage70sq = mwcmc['lage70sq']
lage70sq = np.reshape(lage70sq,(-1,1))

pky50 = mwcmc['pky50']
pky50 = np.reshape(pky50,(-1,1))

munksmk = mwcmc['munksmk'].astype(int)
munksmk = np.reshape(munksmk,(-1,1))

lppd = mwcmc['lppd']
lppd = np.reshape(lppd,(-1,1))

xdgy = mwcmc['xdgy']
xdgy = np.reshape(xdgy,(-1,1))

presur0 = mwcmc['presur0'].astype(int)
presur0 = np.reshape(presur0,(-1,1))


presur1 = mwcmc['presur1'].astype(int)
presur1 = np.reshape(presur1,(-1,1))


presur2 = mwcmc['presur2'].astype(int)
presur2 = np.reshape(presur2,(-1,1))


presur3 = mwcmc['presur3'].astype(int)
presur3 = np.reshape(presur3,(-1,1))


presur4 = mwcmc['presur4'].astype(int)
presur4 = np.reshape(presur4,(-1,1))


presur5 = mwcmc['presur5'].astype(int)
presur5 = np.reshape(presur5,(-1,1))


presur6 = mwcmc['presur6'].astype(int)
presur6 = np.reshape(presur6,(-1,1))


fpusur6 = mwcmc['fpusur6'].astype(int)
fpusur6 = np.reshape(fpusur6,(-1,1))

pudgy = mwcmc['pudgy']
pudgy = np.reshape(pudgy,(-1,1))
mpudgy = male*pudgy
fpudgy = female*pudgy

pytot = 10000 * np.sum(pyr)


print '\n', len(np.unique(udbid)) , 'workers with','%9.0f' % pytot, 'person-years and', '%4.0d' % np.sum(cases),'cases in full cohort'



print "CIM computation model:", models[umod]['title']
# create pandas dataframemodcovs for covariate info
 
stcovs = []
# baseline covariates (LOGL 0)
stcovs.append(np.column_stack((con, female, lage70, lage70sq)))

# smoking covariates Line 1 + LOGL 1
stcovs.append(np.column_stack((pky50, munksmk)))
stcovs.append(np.row_stack((lppd)))
stcovs.append(np.row_stack((xdgy)))

# Pu dose 
if (umod == 2):
    stcovs.append(np.column_stack((mpudgy, fpudgy)))
    pn1 = "Male Pu ERR/Gy"
    pn2 = "Female Pu ERR/Gy"
elif (umod == 0 ):
    stcovs.append(np.row_stack((pudgy)))
    stcovs.append(np.row_stack((female)))
    pn1 = "Male Pu ERR/Gy"
    pn2 = "F:M Pu ERR ratio"
elif (umod == 1 ):
    stcovs.append(np.row_stack((pudgy)))
    stcovs.append(np.row_stack((male)))
    pn1 = "Female Pu ERR/Gy"
    pn2 = "M:F Pu ERR ratio"    

# external dose LINE 2

stcovs.append(np.column_stack((presur0, presur1, presur2, presur3, presur4, presur5, presur6, fpusur6)))


print '\nAnalysis model:', models[modidx]['title'],'with',len(models[modidx]['subterms'])
print "Subterm covariate lists:", len(stcovs)    


mfmodel = models[umod]  # compute subterm fitted values

mfit = modeleval(mfmodel, stcovs, pyr[:, 0])
# compute fitted value vectors
fv = calc_fvalues(mfmodel, mfit['termvals'])
 
showfittedcases(mfmodel,cases,fv,tnames)

mu = np.reshape(fv['fv'],(-1,1))
resid = cases - mu

Q = mkQmat(stcovs, mu, mfit, mfmodel['type'])
colhdr = 'cons,female,lage70,lage70sq,pky50,munksmk,lppd,xdgy,popstpudgy,pufemale,presur0,presur1,presur2,presur3,presur4,presur5,presur6,fpusur6'

# score vector and information
S = Q*resid
parmscores = np.dot(np.transpose(Q),resid)
info = np.dot(S.T,S)

evals = np.linalg.eigvals(info)
print "Information-matrix eigenvalues"
for i in range(len(evals)):
    print '%3d' % i, '%10.5f' % evals[i]
print

# parameter SEs from fit (psv file)
modprmses = np.sqrt(np.diag(mfmodel['prmcov']))

parms = info.shape[0]
prmidx = np.array(range(parms))+1
freeparm = 1 - mfmodel['prmstat'][0:parms]
uparm = prmidx * freeparm
pscore = parmscores[uparm>0]
uparm = uparm[uparm>0]-1

modprmses = modprmses[modprmses>0]

infosub = info[np.ix_(uparm,uparm)]
covmat = np.linalg.inv(infosub)
parmses = np.sqrt(np.diag(covmat))
print "   pno    modse        Qse     Score"
for i in range(len(parmses)):
    print '%5d' % (uparm[i]+1), '%10.5f' % modprmses[i], '%10.5f' %  parmses[i], '%10.5f' % pscore[i][0]
 

prmcov = mfmodel['prmcov']
prmvec = mfmodel['prmvec']
    
doseterm = 3
dp1 = 9
dp2 = 10
    
mfmodel['prmnames'][dp1-1] = pn1
mfmodel['prmnames'][dp2-1] = pn2
 
# compute G-matrix 

nopubase = mu/(1 + mfit['termvals'][doseterm][1])

if umod == 0 :
    dparm = np.exp(prmvec[dp2-1]*female)
elif umod == 1 : #  F + M:F sex ratio
    dparm = np.exp(prmvec[dp2-1]*male) 
elif umod == 2 : # sex-specific ERR
    dparm = np.ones_like(mu) #  * prmvec[dp1]


G = dparm * nopubase
np.savetxt('qmatzz.csv', Q, delimiter=',', header=colhdr, comments='')

np.savetxt('gdiagzz.csv',G, delimiter=',', header='gdiag', comments='')
# compute CIM matrix
cim = mkcimdr(Q, G, intdrep, 1000, prmcov)
np.savetxt('cimzz.csv',cim, delimiter=',', header=colhdr,  comments='')
# Wald and Adjusted bounds
mpuno = dp1-1
mprmest = prmvec[mpuno]
mprmse = np.sqrt(prmcov[mpuno, mpuno])

mmse = np.sqrt(cim[mpuno, mpuno])
mpucim = getCI(mprmest, mprmse, mmse, 0.05)
mpuwald = getCI(mprmest, mprmse, 0, 0.05)
dispcimbnds(pn1,mprmest,mpuwald, mpucim)
modcimprmbnds(mfmodel, cim, level=95)