# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 15:25:55 2017

@author: Dale
"""

import h5py
import os
import numpy as np

trch5dir = 'k:/uralsdosimetry/techa/hdf5'
h5base = 'trc16mc'
h5sub = 'ds1'
h5rz  = 'rz1'


os.chdir(trch5dir)

trch5pg = h5py.File(h5base+h5sub+'.h5','r')
trch5rz = h5py.File("x:/"+h5base+h5rz+'.h5')

trcmembers = trch5pg.keys()

nreps = 1500
npeople = len(trcmembers)
yrfrom = trch5pg.attrs.get('yrfrom')
yrto = trch5pg.attrs.get('yrto')
yrcnt = yrto - yrfrom + 1
orglist = trch5pg.attrs.get('organs')


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

badcnt = 0
badrow = 999999
dzero = np.zeros(66)
badlist = []
# loop over organs
for orgno in range(len(orglist)):  
   organ = orglist[orgno]
   print organ,
   # initialize realization data arrayfor this organ
   rlzdata = np.zeros((npeople,yrcnt,2))        
   # loop over people for this realization
   for rzno in range(0,5):
        rzstr = str(rzno)
        try:
            trch5rz.create_group(rzstr)
        except:
           pass 
            
        rzgrp = trch5rz[rzstr]
        
        if (rzno % 2)==0:
            print  rzno, 
        for people in trch5pg:
            rlzrow = rlzrows[people]
            # print people, rlzrow, rzno, orgno
            rlzdata[rlzrow,:,0] = trch5pg[people + '/int'][rzno,:,orgno]
            try:
                rlzdata[rlzrow,:,1] = trch5pg[people + '/ext'][rzno,:,orgno]
            except:
                if rzno==0:
                    badcnt = badcnt + 1
                    badrow = min(badrow, rlzrow)
                    badlist.append(people)
                    # print "No external dose", people, rlzrow, badcnt, badrow
        
        try:
            rzgrp.create_dataset(organ,data=rlzdata, dtype='f',
                                compression='gzip', compression_opts=9)
        except:
            pass

print "Bad count ", badcnt    
        

'''
rlzdata = np.zeros((npeople,yrcnt,2)) 
pid = '43'        
p43row = rlzrows[pid]  
rzno = 0
orgno= 3        
p43rbmi = trch5pg[pid+'/int'][rzno,:,orgno]
p43rbmx = trch5pg[pid+'/ext'][rzno,:,orgno]


mxrow = -1
for people in trch5pg:
    mxrow = max(mxrow, rlzrows[people])
    print people, rlzrows[people]
    
rlzdata[rlzrow,:,0] = trch5pg[people + '/int'][rzno,:,orgno]
,,,
