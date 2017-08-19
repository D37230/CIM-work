
# coding: utf-8



# load some standard libraries
import os
import sys
import time
import numpy as np   # numerical methods
import h5py  # Hdf5 library
import pandas as pd


# # Define mygithub and hdf5path directories and all else should be fine




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



ungfname = 'lunpyung.csv'
mwcung = np.genfromtxt(modinfodir + ungfname, delimiter=',', names=True)

print "Ungrouped record count:", len(mwcung)



# ## Define subset selection variables
reaux = mwcung['mxplant83'] <= 4
premon2 =(mwcung['tsmcat'] < 3) + 0
hasmon = premon2 == 0
hasmonre= reaux |  ( premon2==0)


print  "Monitored only     ",np.sum(hasmon), 
print  "\nMonitored or Re/Aux", np.sum(hasmonre)

# define index to link nonzero doses to rows in dataset
monindx = np.reshape(range(len(mwcung)),(-1,1))
monrownos  = monindx[hasmon]



# ## Set model form and analysis subset


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



# define index to link nonzero doses to rows in dataset
monindx = np.reshape(range(len(mwcung)),(-1,1))
monrownos  = monindx[hasmon]

mwcungsub = mwcung[hasmon]
# extract subset of rows from ungrouped data file for
# postmonitored years for monitored workers


# get interplatoino info: bounding columns in hdf5 dose array, fractional year
lag = 5
imdinfo = setdoseintinfo(ipindx, mwcung, 'udbid', 'yr', lag, ifromyr, itoyr)
imdinfosub = setdoseintinfo(ipindx, mwcung[hasmon], 'udbid', 'yr', lag, ifromyr, itoyr)
idbar = doseinterp(irlzdb, 'sumstat', 'doses/am', imdinfo, scale=dscale, dlayer=lngidx)
idbarsub = doseinterp(irlzdb, 'sumstat', 'doses/am', imdinfosub, scale=dscale, dlayer=lngidx)

# ## Define model covariates

mwcungsub = mwcung

udbid = mwcungsub['udbid']


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
year = np.reshape(mwcungsub['yr'],(-1,1))

mxplant83 = np.reshape(mwcungsub['mxplant83'],(-1,1))
reaux = mxplant83 < 5 

pytot = 10000 * np.sum(pyr)
pymon = 10000 * np.sum(pyr[hasmon])

print '\n', len(np.unique(udbid)) , 'workers with','%9.0f' % pytot, 'person-years and', '%4.0d' % np.sum(cases),'cases in full cohort'
 
print  len(np.unique(udbid[hasmon])) , 'workers with', '%9.0f' % pymon, 'person-years and', '%4.0d' % np.sum(cases[hasmon]),'cases for post-monitored workers'


# create pandas dataframemodcovs for covariate info
 
ivlist = ['udbid', 'male', 'female', 'cases', 'munksmk', 'presur0', 'presur1', 'presur2', 'presur3', 'presur4', 'presur5', 'presur6', 'fpusur6'  ]

for var in ivlist:
    v =eval(var).astype(int)    
    var = v.astype(int)
    
dfmodcovs = pd.DataFrame({'rowno': monindx.ravel()})
    
#female = female.astype(int)
#cases = cases.astype(int)



for var in ['udbid', 'cases', 'pyr', 'year', 'male', 'female', 'lage70' , 'lage70sq', 'pky50', 'munksmk', 'lppd', 'xdgy', 'presur0', 'presur1', 'presur2', 'presur3', 'presur4', 'presur5', 'presur6', 'fpusur6', 'pudgy']:
    v = eval(var)
    v = v.ravel()
    print var, v.shape 
    dfmodcovs[var] = pd.Series(v,index=dfmodcovs.index) 

dfmodcovs.keys()

dfmodcovs.to_csv(modinfodir+'modcovs.csv',na_rep='.')

# get inerpolated doses for a set of reps

irlzbase = irlzdb['realizations']
ihasmon = hasmon.astype(int)
ihasmon = np.reshape(ihasmon,(-1,1))

udbid = np.reshape(udbid,(-1,1))
udbid = udbid.astype(int)
dfreps = pd.DataFrame({'rowno': monrownos.ravel(), 'udbid': udbid[ihasmon==1].ravel(), 'year': year[ihasmon==1].ravel(), 'pumean': idbar[ihasmon==1].ravel()})
for  rep in range(1000):
    rstr = str(rep)
    rname = 'rep'+rstr
    repdose = doseinterp(irlzbase, rstr, 'doses', imdinfosub, scale=dscale, dlayer=lngidx)
    repdose = repdose.ravel()
    dfreps[rname] = pd.Series(repdose,index=dfreps.index)
    if ((rep+1) % 100) == 0:
        print "rep no", rep
    

dfreps.to_csv(modinfodir+'dfreps.csv',na_rep='.')