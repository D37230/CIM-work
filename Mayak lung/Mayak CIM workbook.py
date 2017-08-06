
# coding: utf-8

# In[1]~

# load some standard libraries
import os
import sys
import time
import numpy as np   # numerical methods
import h5py  # Hdf5 library


# # Define mygithub and hdf5path directories and all else should be fine

# In[2]~


mygithub = 'C:/Users/Dale/Documents/GitHub/'
projdir = mygithub + 'CIM-work/Mayak lung/'
hdf5path = 'k:/UralsDosimetry/'
modinfodir = hdf5path + 'mayak/mwds16lung/'
rlzdir = hdf5path + 'mayak/mwds2016/'
os.chdir(projdir)
sys.path.append('../pythonutils')

# import utility functions for this applicaiton

from psvdict2 import modeleval, getmodinfo, calc_fvalues, showfittedcases
from doseintfun import setdoseintinfo, doseinterp
from cimfuncs import mkcim, getCI, dispcimbnds, mkQmat, modcimprmbnds



# # Set up hdf5 file

# In[3]~


irlzfn = 'res_ad13(final).h5'
irlzdb = h5py.File(rlzdir + irlzfn, 'r')

# associative array for access to internal dose person-specific info
ipeople_list = irlzdb['sumstat/idlist'][:]  # list of ids in row order
ipeople = len(ipeople_list)
ipindx = {}
for i, prsn in enumerate(ipeople_list):
    ipindx[str(prsn)] = int(i)

# associative array for access to internal dose oragns
iorgs_list = irlzdb['sumstat/orglist'][:]  # list of organs with doses
iorgs = len(iorgs_list)
iorgindx = {}
for i, org in enumerate(iorgs_list):
    iorgindx[org] = int(i)
# earliest internal exposure year (column 1)
ifromyr = irlzdb['sumstat'].attrs['fromyr']
# last exposure year
itoyr = irlzdb['sumstat'].attrs['toyr']
inyrs = itoyr - (ifromyr) + 1

# get dose units and set dose scvalar factor (models are fit using dose in Gy)
dsunits = irlzdb['sumstat/doses'].attrs['units']
dscale = 1
if dsunits == 'mGy':
    dscale = 1000
elif dsunits == 'cGy':
    dscale == 100
elif dsunits == "Gy":
    dscale = 1

# define lung dose layer
lngidx = iorgindx['LUNG']


# In[4]~


# associative array for access to internal dose person-specific info
ipeople_list = irlzdb['sumstat/idlist'][:]  # list of ids in row order
ipeople = len(ipeople_list)
ipindx = {}
for i, prsn in enumerate(ipeople_list):
    ipindx[str(prsn)] = int(i)

# associative array for access to internal dose oragns
iorgs_list = irlzdb['sumstat/orglist'][:]  # list of organs with doses
iorgs = len(iorgs_list)
iorgindx = {}
for i, org in enumerate(iorgs_list):
    iorgindx[org] = int(i)
# earliest internal exposure year (column 1)
ifromyr = irlzdb['sumstat'].attrs['fromyr']
# last exposure year
itoyr = irlzdb['sumstat'].attrs['toyr']
inyrs = itoyr - (ifromyr) + 1

# get dose units and set dose scvalar factor (models are fit using dose in Gy)
dsunits = irlzdb['sumstat/doses'].attrs['units']
dscale = 1
if dsunits == 'mGy':
    dscale = 1000
elif dsunits == 'cGy':
    dscale == 100
elif dsunits == "Gy":
    dscale = 1

# define lung dose layer
lngidx = iorgindx['LUNG']


# # Open ungrouped data file from modinfodir

# In[5]~

ungfname = 'lunpyung.csv'
mwcung = np.genfromtxt(modinfodir + ungfname, delimiter=',', names=True)

print "Ungrouped record count:", len(mwcung)


# ## Define subset selection variables 

# In[6]~

reaux = mwcung['mxplant83'] <= 4
premon2 =(mwcung['tsmcat'] < 3) + 0
hasmon = premon2 == 0
hasmonre= reaux |  ( premon2==0)

print  "Monitored only     ",np.sum(hasmon), 
print  "\nMonitored or Re/Aux", np.sum(hasmonre)


# ## Set model form and analysis subset

# In[7]~

'''
    umod = 0 | 3   # male + F:M sex ratio model
    umod = 1 | 4   # female + M"F sex ratio model
    umod = 2 | 5   # male +_ female dose effect model
    
    msub = 0 # all data with Pu Surrogate
    msub = 1 # Monitored + Re/Aux
    msub = 2 # Montored workers only

'''
msub = 0
umod = 0
psvname = 'mwds16pualtmods.psv'

tnames = ["Baseline", "Smoking", "External dose", "Pu Surrogate" ,"Pu dose"]
if msub == 1:
    mwcungsub = mwcung[hasmonre]
    modidx = umod
elif msub == 2:
    mwcungsub = mwcung[hasmon]
    umod = umod
    modidx = umod+3
else:
    mwcungsub = mwcung
    psvname = 'mwds16pu.psv' 
    modidx = umod
    tnames = ["Baseline", "Smoking", "External dose", "Pu Surrogate" ,"Pu dose"]
    
print len(mwcungsub), 'records in anlaysis subset'


# ## Read PSV file

# In[8]~

psavef = open(modinfodir + psvname, 'r')
models = getmodinfo(psavef)
for mno in range(len(models)):
    print "model ",mno,models[mno]['title']
psavef.close()

print "CIM computation model:", models[umod]['title']


# ## Get basic info for internal dose interpolation and interpolate mean doses

# In[9]~

lag = 5
imdinfo = setdoseintinfo(ipindx, mwcungsub, 'udbid', 'yr', lag, ifromyr, itoyr)
idbar = doseinterp(irlzdb, 'sumstat', 'doses/am', imdinfo, scale=dscale, dlayer=lngidx)


# ## Define model covariates

# In[10]~


male = np.reshape(((mwcungsub['sex'] == 0) + 0),(-1,1))
female = 1 - male
age = np.reshape(mwcungsub['age'],(-1,1))
lage70 = np.log(age /70)
lage70sq = lage70*lage70
unksmk = np.reshape((mwcungsub['smkcat'] == 1) + 0,(-1,1))
munksmk = male * unksmk
smoker = np.reshape(((mwcungsub['smkcat'] == 3) + 0),(-1,1))
neversmk = np.reshape(((mwcungsub['smkcat'] == 2) + 0),(-1,1))
smkppd =np.reshape(mwcungsub['ppd'],(-1,1)) * (1 - unksmk)
smkdur =  (age - np.reshape(mwcungsub['smkstart'],(-1,1)))* smoker
smkdur = smkdur * (smkdur > 0)
packyrs = smkppd * smkdur * smoker
pky50 = packyrs/50
lppd = np.log(smkppd * smoker + unksmk + neversmk)

 
premon2 = (np.reshape(mwcungsub['tsmcat'],(-1,1)) < 2) + 0
postmon2 = 1 - premon2
pudgy = np.reshape((mwcungsub['lag5id']/1000),(-1,1)) 

xdgy =np.reshape(( mwcungsub['lag5ed']/1000),(-1,1)) 
postpudgy = pudgy * postmon2
mpudgy = male * postpudgy
fpudgy = female * postpudgy

presur = premon2 * np.reshape(mwcungsub['pusur'],(-1,1))
presur0 = presur == 0
presur1 = presur == 1
presur2 = presur == 2
presur3 = presur == 3
presur4 = presur == 4
presur5 = presur == 5
presur6 = presur == 6
fpusur6 = female * presur6

pyr = np.reshape((np.array(mwcungsub['pyr']/10000.)),(-1,1))
cases = np.reshape(mwcungsub['lung'],(-1,1))

mxplant83 = np.reshape(mwcungsub['mxplant83'],(-1,1))
reaux = mxplant83 < 5 

print "Total cases: ", '%6.0d' % np.sum(cases)


# ## Define model-specific subterm covariate vector list

# In[11]~

stcovs = []
# baseline covariates (LOGL 0)
stcovs.append(np.column_stack((male, female, lage70, lage70sq)))

# smoking covariates Line 1 + LOGL 1
stcovs.append(np.column_stack((pky50, munksmk)))
stcovs.append(np.row_stack((lppd)))

# external dose LINE 2
stcovs.append(np.row_stack((xdgy)))

if msub==0:
    # pusur (preemonitoring)
    stcovs.append(np.column_stack((presur0, presur1, presur2, presur3, presur4, presur5, presur6, fpusur6)))

# Pu dose 
if (umod == 2):
    stcovs.append(np.column_stack((mpudgy, fpudgy)))
    
    
elif (umod == 0 ):
    stcovs.append(np.row_stack((pudgy)))
    stcovs.append(np.row_stack((female)))
elif (umod == 1 ):
    stcovs.append(np.row_stack((pudgy)))
    stcovs.append(np.row_stack((male)))

print '\nAnalysis model:', models[modidx]['title'],'with',len(models[modidx]['subterms'])
print "Subterm covariate lists:", len(stcovs)    



# ## Fitted value computaitons for this model in these data

# In[12]~


mfmodel = models[umod]  # compute subterm fitted values

mfit = modeleval(mfmodel, stcovs, pyr[:, 0])
# compute fitted value vectors
fv = calc_fvalues(mfmodel, mfit['termvals'])
 
showfittedcases(mfmodel,cases,fv,tnames)

mu = np.reshape(fv['fv'],(-1,1))
resid = cases - mu



# ## Computation of Q matrix
# $$ Q = \frac{ d \mu / d \beta}{\mu }$$

# In[13]~


Q = mkQmat(stcovs, mu, mfit, mfmodel['type'])


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
    



# ## G matrix for various Pu models  $$   \frac{d \mu}{d x}$$

# In[14]~

prmcov = mfmodel['prmcov']
prmvec = mfmodel['prmvec']

haspudose =np.reshape(( reaux | (premon2==0)),(-1,1))

doseterm = 3
dp1 = 9
dp2 = 10
if msub == 0:
    doseterm = 4
    dp1 = 17
    dp2 = 18
    
if umod == 0:
    mfmodel['prmnames'][dp1-1] = 'Male Pu ERR'
    mfmodel['prmnames'][dp2-1] = 'F:M Pu ERR sex ratio'
elif umod == 1:
    mfmodel['prmnames'][dp1-1] = 'Female Pu ERR'
    mfmodel['prmnames'][dp2-1] = 'M:F Pu ERR sex ratio'
elif umod == 2:
    mfmodel['prmnames'][dp1-1] = 'Male Pu ERR'
    mfmodel['prmnames'][dp2-1] = 'Female Pu ERR'   

nopubase = mu/(1 + mfit['termvals'][doseterm][1])

if umod == 0 :
    dparm = np.exp(prmvec[dp2-1]*female)
elif umod == 1 : #  F + M:F sex ratio
    dparm = np.exp(prmvec[dp2-1]*male) 
elif umod == 2 : # sex-specific ERR
    dparm = np.ones_like(mu) #  * prmvec[dp1]
    # dparm[female] = prmvec[dp2] 
    # dparm = dparm * haspudose


G = dparm * nopubase


# ## CIM setup

# In[15]~

irlzbase = irlzdb['realizations']
rlzcnt = len(irlzbase.keys())

prmcov = mfmodel['prmcov']
prmvec = mfmodel['prmvec']

irlzbase = irlzdb['realizations']
rlzcnt = len(irlzbase.keys())

prmcov = mfmodel['prmcov']
prmvec = mfmodel['prmvec']


cim = mkcim(Q, G, irlzbase, rlzcnt, 'doses', imdinfo, prmcov, dslayer=lngidx, scale = dscale)


# ## Pu ERR CIM

# In[16]~


mpuno = dp1-1
mprmest = prmvec[mpuno]
mprmse = np.sqrt(prmcov[mpuno, mpuno])

mmse = np.sqrt(cim[mpuno, mpuno])
mpucim = getCI(mprmest, mprmse, mmse, 0.05)
mpuwald = getCI(mprmest, mprmse, 0, 0.05)
dispcimbnds('Female Pu dose response ERR:',mprmest,mpuwald, mpucim)


# In[18]~

modcimprmbnds(mfmodel, cim, level=95)
modcimprmbnds(mfmodel, cim, level=90)
modcimprmbnds(mfmodel, cim, level=68)
modcimprmbnds(mfmodel, cim, level=99)


# In[ ]~



