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
    
eurtmembers = eurth5pg.keys()

nreps = 1500
npeople = len(eurtmembers)
yrfrom = eurth5pg.attrs.get('yrfrom')
yrto = eurth5pg.attrs.get('yrto')
yrcnt = yrto - yrfrom + 1
orglist = eurth5pg.attrs.get('organs')
dlayers = eurth5pg.attrs.get('dlayers')



# create sorted list of id's for realization data set
idlist = []
for chid in eurtmembers:
    idlist.append(int(chid))

idlist = np.sort(idlist)

# create dictionary for rows in realiztion datasets
rlzrows= {}
rowno = 0
for id in idlist:
    rlzrows[str(id)] =  rowno
    rowno = rowno + 1
 
orgindx = {}
for i in range(len(orglist)):
    orgindx[orglist[i]] = i


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