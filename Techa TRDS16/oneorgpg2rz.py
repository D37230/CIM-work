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

trch5pg = h5py.File(h5base+h5sub+'fix.h5','r')
trch5rz = h5py.File("x:/"+h5base+h5rz+'test2.h5')
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
yrfrom = 1950
yrto = 2015
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


# set uyp for making realizaiton data sets
    

badcnt = 0
badrow = 999999
badlist = []
# loop over organs
organ = 'sto'
orgno = 2
try:
   ssorg = ssgrp.create_group(organ)
except:
   ssorg = ssgrp[organ]
print '\n',organ,
chunksz = 40
nchunks = nreps/chunksz
# initialize realization data arrayfor this organ
rlzdata = np.zeros((chunksz,npeople,yrcnt,2))
rzam = np.zeros((npeople,yrcnt,2))         
# loop over people for this realization
wgt=0
for chunk in [0]:
        rznos = range(chunk*chunksz,(chunk+1)*chunksz)
        for rzno in rznos:
            rzstr = str(rzno)
            rzngrp = rzgrp.create_group(rzstr)
        
            
        i = 0
        print '\nChunk ',chunk,':',
        for people in trch5pg.iteritems():
            i += 1
            if i % 5000 == 0:
                print i,
            rlzrow = rlzrows[people[0]]
            # print people, rlzrow, rzno, orgno
            rzdi = people[1]['int'][rznos,:,orgno] 
            rlzdata[range(0,chunksz),rlzrow,:,0] = rzdi
            rzdx = people[1]['ext'][rznos,:,orgno]
            rlzdata[range(0,chunksz),rlzrow,:,1] = rzdx

        for i in range(0,chunksz):
            rzngrp = rzgrp[str(rznos[i])]
            rzngrp.create_dataset(organ,data=rlzdata[i,:,:,:], dtype='f',
                                  compression='gzip', compression_opts=9)    
'''    

        wgt += 1
        rzam = rzam + (rlzdata - rzam)/wgt
        
ssorg.create_dataset('am',data=rzam, compression='gzip',compression_opts=9)
'''        

   
        

trch5rz.close()
