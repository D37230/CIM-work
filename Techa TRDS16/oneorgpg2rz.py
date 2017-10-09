# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 15:25:55 2017

@author: Dale
"""

import h5py
import os
import sys
import numpy as np

ghdir = 'C:/Users/Dale/Documents/GitHub/CIM-Work/'
cimworkdir = ghdir + 'Techa TRDS16'
sys.path.append(cimworkdir)

import hdf5utils as trch5


trch5dir = 'k:/uralsdosimetry/techa/hdf5/'
h5base = 'trc16mc'
h5sub = 'ds1'
h5rz  = 'trz1'

pgname = trch5dir +h5base+h5sub+'fix.h5'
rzname = 'x:/' + h5base +h5rz+'.h5'
os.chdir(trch5dir)

trch5pg = h5py.File(pgname,'r')
trch5rz = h5py.File(rzname)

# copy attributes to root in realization file
for atr in trch5pg.attrs.keys():
    trch5rz.attrs[atr] =  trch5pg.attrs[atr]
    print atr, trch5rz.attrs[atr]
    
    
try:
    rzgrp = trch5rz.create_group('realizations')
except:
     rzgrp = trch5rz['realizations']   
try:
    ssgrp = trch5rz.create_group('sumstat')
except:
    ssgrp = trch5rz['sumstat']
   
#  get file attributes
trcmembers = trch5pg.keys()
npeople = len(trcmembers)

nreps = trch5pg.attrs['nreps']

yrfrom = trch5pg.attrs['yrfrom']
yrto = trch5pg.attrs['yrto']
yrcnt = yrto - yrfrom + 1
orglist = trch5pg.attrs['organs']


# create sorted list of id's for realization data set
idlist = []
for chid in trcmembers:
    idlist.append(int(chid))

idlist = np.sort(idlist)

# create dictionary for rows in realiztion datasets
rlzrows= {}
rowno = 0
for id in idlist:
    rlzrows[str(id)] =  rowno
    rowno = rowno + 1


# set uyp for making realizaiton data sets
    

# loop over organs


for orgno in range(len(orglist)):
    organ = orglist[orgno]
    print orgno, organ 
    trch5.pgtorlz(trch5pg, rzgrp, npeople, nreps, yrcnt, orgno, organ, rlzrows)   
        
trch5pg.close()


rzgrp = trch5rz['realizations']
for orgno in range(len(orglist)):
    rzam = np.zeros((npeople,yrcnt,2))
    rzse = rzam
    rno = 0
    organ = orglist[orgno]

    print '\nComputing mean', organ, 'dose:', 
    for rlz in rzgrp.iteritems():
        rno += 1
        rzdata = rlz[1][organ][:]
        delta = (rzdata - rzam)/rno
        rzam = rzam + delta
        rzse = rzse + (rzdata-rzam)*delta
        if (rno % 100)==0:
            print rno,
    rzse = rzse/max(1,(rno-1))
 
    try:
        ssorg = trch5rz.create_group('sumstat/'+organ)
    except:
        ssorg = trch5rz['sumstat/'+organ]
 

    ssorg.create_dataset('am',data=rzam, dtype='f',
                          compression='gzip', compression_opts=9) 
    # ssorg.create_dataset('/se',data=rzse, dtype='f',
    #                      compression='gzip', compression_opts=9) 

trch5rz.close() 