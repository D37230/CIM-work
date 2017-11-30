# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 15:25:55 2017

@author: Dale
"""

import h5py
import os
import numpy as np

eurth5dir = 'h:/uraldosimetry/techa/hdf5'
h5base = 'eurt16mc'
h5sub = ''
h5rz  = 'rz1'


os.chdir(eurth5dir)

eurth5pg = h5py.File(h5base+h5sub+'.h5','r')
eurth5rz = h5py.File(h5base+h5rz+'.h5')
try:
    rzgrp = eurth5rz.create_group('realizations')
except:
     rzgrp = eurth5rz['realizations']   
try:
    ssgrp = eurth5rz.create_group('sumstat')
except:
    ssgrp = eurth5rz['sumstat']
    
eurtmembers = eurth5pg.keys()

nreps = 1500
npeople = len(eurtmembers)
yrfrom = eurth5pg.attrs.get('yrfrom')
yrto = eurth5pg.attrs.get('yrto')
yrcnt = yrto - yrfrom + 1
orglist = eurth5pg.attrs.get('organs')
dlayers = eurth5pg.attrs.get('dlayers')

eurth5rz.attrs.modify('yrfrom',yrfrom)
eurth5rz.attrs.modify('yrto',yrto)
eurth5rz.attrs.modify('dlayers',dlayers)


# create sorted list of id's for realization data set
idlist = []
for chid in eurtmembers:
    idlist.append(int(chid))

idlist = np.sort(idlist)

try:
    ssgrp.create_dataset('idlist',data=idlist, dtype='l',
                          compression='gzip', compression_opts=9)
except:
    pass

# create dictionary for rows in realiztion datasets
rlzrows= {}
rowno = 0
for id in idlist:
    rlzrows[str(id)] =  rowno
    rowno = rowno + 1
    
    
'''
# count leading 0 doses    
dzeros = {}


npeeps=0
for people in eurth5pg :
    npeeps = npeeps+1
    if npeeps%1000 == 0:
        print npeeps,
    r0 = eurth5pg[people + '/int'][0,:,0]
    lzcnt  = 0    
    for i in range(0,yrcnt):
        if r0[i] == 0.:
            lzcnt = lzcnt+ 1
        else:
            dzeros[people]=lzcnt
'''

# set uyp for making realizaiton data sets
    

badcnt = 0
badrow = 999999
badlist = []

rz1orgnos = [3]
rz1orgs = orglist[[3]]
# loop over organs
for orgno in rz1orgnos: 
   organ = orglist[orgno]
   try:
        ssorg = ssgrp.create_group(organ)
   except:
        ssorg = ssgrp[organ]
   print '\n',organ,
   # initialize realization data arrayfor this organ
   rlzdata = np.zeros((npeople,yrcnt,2))
   # loop over people for this realization
   wgt=0
   for rzno in range(0,nreps):
        rzstr = str(rzno)
        try:
           rzngrp = rzgrp.create_group(rzstr)
        except:
           rzngrp = rzgrp[rzstr]
        
        if (rzno % 100)==0:
            print  rzno, 
        for people in eurth5pg:
            pstr = str(people)
            rlzrow = rlzrows[people]
            # print people, rlzrow, rzno, orgno
            # lz = dzeros[people]
            rzdi = eurth5pg[pstr + '/int'][rzno,:,orgno] 
            # if lz>0:
            #    rzdi[range(0,lz)] = 0
            rlzdata[rlzrow,:,0] = rzdi
            try:
                rzdx = eurth5pg[pstr + '/ext'][rzno,:,orgno]
                # if lz>0:
                #   rzdx[range(0,lz)] = 0
                rlzdata[rlzrow,:,1] = rzdx
            except:
                if rzno==0:
                    badcnt = badcnt + 1
                    badrow = min(badrow, rlzrow)
                    badlist.append(people)
                    # print "No external dose", people, rlzrow, badcnt, badrow
        
        try:
            rzngrp.create_dataset(organ,data=rlzdata, dtype='f',
                                compression='gzip', compression_opts=9)
        except:
            pass
          

print "\nBad count ", badcnt    
        
eurth5rz.close()