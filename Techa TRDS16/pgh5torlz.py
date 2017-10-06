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
try:
    rzgrp = trch5rz.create_group('realizations')
except:
     rzgrp = trch5rz['realizations']   
try:
    ssgrp = trch5rz.create_group('sumstat')
except:
    ssgrp = trch5rz['sumstat']
    
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
# count leading 0 doses    
dzeros = {}
npeeps=0
for people in trch5pg :
    npeeps = npeeps+1
    if npeeps%1000 == 0:
        print npeeps,
    r0 = trch5pg[people + '/int'][0,:,0]
    lzcnt  = 0    
    for i in range(0,yrcnt):
        if r0[i] == 0.:
            lzcnt = lzcnt+ 1
        else:
            dzeros[people]=lzcnt
            break
print

# set uyp for making realizaiton data sets
    

badcnt = 0
badrow = 999999
badlist = []
# loop over organs
for orgno in [0,2]:  
   organ = orglist[orgno]
   try:
        ssorg = ssgrp.create_group(organ)
   except:
        ssorg = ssgrp[organ]
   print '\n',organ,
   # initialize realization data arrayfor this organ
   rlzdata = np.zeros((npeople,yrcnt,2))
   rzam = np.zeros((npeople,yrcnt,2))         
   # loop over people for this realization
   wgt=0
   for rzno in range(0,1500):
        rzstr = str(rzno)
        try:
           rzngrp = rzgrp.create_group(rzstr)
        except:
           rzngrp = trch5rz[rzstr]
        
        if (rzno % 100)==0:
            print  rzno, 
        for people in trch5pg:
            rlzrow = rlzrows[people]
            # print people, rlzrow, rzno, orgno
            lz = dzeros[people]
            rzdi = trch5pg[people + '/int'][rzno,:,orgno] 
            if lz>0:
                rzdi[range(0,lz)] = 0
            rlzdata[rlzrow,:,0] = rzdi
            try:
                rzdx = trch5pg[people + '/ext'][rzno,:,orgno]
                if lz>0:
                    rzdx[range(0,lz)] = 0
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
            wgt += 1
            rzam = rzam + (rlzdata - rzam)/wgt
        except:
            pass
   try:
       ssorg.create_dataset('am',data=rzam, compression='gzip',
                             compression_opts=9)
   except:
        print '\n could not create ',organ,' am sumstat dataset'
        

print "\nBad count ", badcnt    
        


'''

pid = '43'        
p43row = rlzrows[pid]  
rzno = 0
orgno= 3        
p43rbmi = trch5pg[pid+'/int'][rzno,:,orgno]
p43rbmx = trch5pg[pid+'/ext'][rzno,:,orgno]

r0rbm = trch5rz['0/rbm'][:]
r1rbm = trch5rz['1/rbm'][:]


mxrow = -1
for people in trch5pg:
    mxrow = max(mxrow, rlzrows[people])
    print people, rlzrows[people]
    
rlzdata[rlzrow,:,0] = trch5pg[people + '/int'][rzno,:,orgno]
'''