# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 10:48:27 2017

@author: prest
"""

import os
import sys
import numpy as np   # numerical methods
import h5py  # Hdf5 library
import pandas as pd

dsrlzindir = "h:/Uraldosimetry/techa/napier/"
dsrlzh5dir = "w:/russia/urcrmcurrent/TRDS16MC/"

eurth5fn = 'eurtds16mc.h5'
trch5fn = 'trcds16mc.h5'



eurticsv = "Development.EIN"
eurtxcsv = "Development.EOU"
trcicsv  = "Development.TIN"
trcxcsv  = "Development.TOU"

os.chdir(dsrlzh5dir)

eurth5pg = h5py.File(dsrlzh5dir+eurth5fn)

eurtihdr = pd.read_csv(dsrlzindir+eurticsv, nrows=0)
eurtihdr = eurtihdr.rename(columns = {'ubid':'sysnum'})

# get dose column names and sert dose names
dscols =  eurtihdr.columns.tolist()[3:] 
dsnames  =  [dn.replace('ein','') for dn in dscols] 
dscnt = len(dsnames)

# read all rows except last one
eurtidat = pd.read_csv(dsrlzindir+eurticsv, header=None, skiprows=[0], 
                       sep='\s+', names = eurtihdr)[:-1]
eurtxdat = pd.read_csv(dsrlzindir+eurtxcsv, header=None, skiprows=[0], 
                       sep='\s+', names = eurtihdr)[:-1]

people = pd.unique(eurtidat['sysnum'])
peoplecnt = len(people)
yrcnt = len(pd.unique(eurtidat['year']))
yrfrom = np.min(eurtidat['year'])
yrto = np.max(eurtidat['year'])
rlzcnt = len(pd.unique(eurtidat['real']))
rlzmax = np.max(eurtidat['real'])

print '\n', peoplecnt,'people with', 
print  yrcnt, 'years from ', yrfrom, 
print 'to',yrto,'with', dscnt, 'organ doses'

if rlzcnt == rlzmax:
    print rlzcnt, 'realizations per person-year\n'
else: 
    print '\n\nNumber of unique realizations (', rlzcnt, 
    print ') and maximum realization number',rlzmax,'are inconsistent\n'
    
# sort data frame by sysnum, realization, year
eurtidat = eurtidat.sort_values(by=['sysnum', 'real', 'year'])
eurtxdat = eurtxdat.sort_values(by=['sysnum', 'real', 'year'])
# get doses for one person
def addtoh5pg(h5pg, p, pu, px, dsdat, dstype, dscols):
    '''
    dosestoh5pg - save doses in csv file to HDF5 file with person groups
    
    Arguments:
            h5pg
                hdf5file for doses
            
            p
                person identifier
                
            dsdat
                pandas data fream with dose realization data
                (one row per person-year) with columns sysnum, real, year, and doses
            
            dstype
                'int' for internal or 'ext' for external
            
            dscols
                list of dose column names
    '''
    doses = dsdat.loc[dsdat['sysnum']==p, 'real':dscols[-1]]
    # fill 0 dose with previous years value
    doses = doses.apply(lambda x: x.replace(to_replace=0, method='ffill'), axis=0)
    dsarr = doses.as_matrix(columns=dscols)
    dsarr = dsarr.reshape((rlzcnt,yrcnt,-1),order='C')
    dsmn = np.mean(dsarr,axis=0)
    pgroup = pstr + '/' + dstype
    try:
        h5pg[pgroup][:]
        px += dstype + ' '
    except:
        pu += dstype + ' ' 
        h5pg.create_dataset(pgroup,data=dsarr,dtype='f',compression='gzip',compression_opts=9)
    return pu, px, dsmn

people = pd.unique(eurtidat['sysnum'])
try:
    h5pg = eurth5pg.create_group('people')
except:
    h5pg = eurth5pg['people']
    
for p in people:
    pstr = str(p)
    print p,
    px = ''
    pu = ''
    pu, px, idsmn = addtoh5pg(h5pg, p, pu, px, eurtidat, 'int', dscols)
    pu, px, xdsmn = addtoh5pg(h5pg, p, pu, px, eurtxdat, 'ext', dscols)
    dsmn = np.dstack((idsmn, xdsmn))
    pgroup = pstr + '/dsmn'
    try:
        h5pg[pgroup][:]
        px += 'dsmn'
    except: 
        h5pg.create_dataset(pgroup,data=dsmn,dtype='f',compression='gzip',compression_opts=9)
        pu += 'dsmn'
    if pu:
        print pu,'added',
    if px:
        print ':',px,'exist',
    print
        
        
eurth5pg.attrs['yrfrom']=yrfrom    
eurth5pg.attrs['yrto']=yrto
eurth5pg.attrs['organs'] = dsnames
eurth5pg.attrs['units'] = 'mGy'
eurth5pg.attrs['dlayers'] = ['int', 'ext']
   
eurth5pg.close()
