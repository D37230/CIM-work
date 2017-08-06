# -*- coding: utf-8 -*-
"""
Created on Tue Aug 01 10:58:47 2017

@author: Dale
"""

# In[1]:

# load some standard libraries
import os
import sys
import time
import numpy as np   # numerical methods
import h5py  # Hdf5 library


# In[2]:

# move to working directory
# redine mygithub path and all else should be fine

mygithub = 'C:/Users/Dale/Documents/GitHub/'
projdir = mygithub + 'CIM work/Mayak lung/'
hdf5path = 'k:/UralsDosimetry/'
os.chdir(projdir)
sys.path.append('../pythonutils')

# functions for this applicaiton

from psvdict2 import modeleval, getmodinfo, calc_fvalues
from doseintfun import setdoseintinfo, doseinterp
from cimfuncs import mkcim, getCI, dispcimbnds, mkQmat



# In[3]:

# directories for model info and dose realization HDF5 files
rlzdir = hdf5path + 'mayak/mwds2016/'
# set up hdf5 file

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

# define lung dose layer
lngidx = iorgindx['LUNG']

# In[4]:

yrnames = ['udbid']
for yr in range(itoyr-ifromyr+1):
    yrnames.append('yr'+str(ifromyr + yr))

yrnames= np.array(yrnames)
yrnames= np.reshape(yrnames,(1,-1))

mnfile = open(rlzdir+'mwds16lunam.csv','w')
np.savetxt(mnfile,yrnames,delimiter=',',fmt='%s') 

lunmean = irlzdb['sumstat/doses/am'][:,:,lngidx]
lunmean = np.column_stack((ipeople_list, lunmean))
np.savetxt(mnfile,lunmean,delimiter=',') 

mnfile.close()

# In[15]:

yrnames = ['rlzno', 'udbid']
for yr in range(itoyr-ifromyr+1):
    yrnames.append('yr'+str(ifromyr + yr))

yrnames= np.array(yrnames)
yrnames= np.reshape(yrnames,(1,-1))

rlzfile = open(rlzdir+'mwds16lunrlz.csv','w')
rlzno = np.zeros_like(lunmean[:,0])
rlzno = np.reshape(rlzno,(-1,1))
np.savetxt(rlzfile,yrnames,delimiter=',',fmt='%s') 

rlzno
for i in range(1000):

    rlzstr = 'realizations/' + str(i) +'/doses'
    print np.mean(rlzno), rlzstr, rlzno.shape
    lunrlz = irlzdb[rlzstr][:,:,lngidx]
    lunrlzaug = np.column_stack((ipeople_list, lunrlz))
    lunrlzaug = np.column_stack((rlzno,lunrlzaug))
    np.savetxt(rlzfile,lunrlzaug,delimiter=',')
    rlzno = rlzno + 1
    
rlzfile.close()
