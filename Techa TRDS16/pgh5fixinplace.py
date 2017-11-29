# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 15:25:55 2017

@author: Dale
"""

import h5py
import os
import numpy as np

'''
  Fix zero dose problems in TRC person group files
  a)  initial zero doses for realizations after the first one were incorrectly filled
      with value of total cumulative dose for previous realization
  b)  datasets for 260 TRC members with non-zero internal doses but 0 external doses 
      were not written to the hdf5 file
  c)  individual with 0 internal and external doses was not included in file (sysnum 132119)

'''

trch5dir = 'k:/uralsdosimetry/techa/hdf5'
trch5fxdir = 'f:/techa_hdf5'
h5base = 'trc16mc'
h5sub = 'ds3'

os.chdir(trch5fxdir)


trch5pg = h5py.File(h5base+h5sub+'.h5','r+')
trch5fx = h5py.File(trch5fxdir + '/' + h5base+h5sub+'fix.h5')

pgname = trch5dir + '/' + h5base + h5sub + '.h5'
fxname = trch5fxdir + '/' + h5base + h5sub + 'fix.h5'
# add group for person with 0 internal and external MC dose esitmates
try:
    trch5pg['132119']
    print '132119 already exists'
    print trch5pg['132119'].keys()

except:
    trch5pg.create_group('132119')
    
    
trcmembers = trch5pg.keys()
trch5pg.attrs.modify('nreps',1500)
trch5pg.attrs.modify('npeople',len(trcmembers))
# reopen as read only
trch5pg.close()
trch5pg = h5py.File(trch5dir + '/'+ h5base+h5sub+'.h5','r')

# get attributes from old file
trcmembers = trch5pg.keys()
nreps = trch5pg.attrs.get('nreps')
npeople = trch5pg.attrs.get('npeople')
yrfrom = trch5pg.attrs.get('yrfrom')
yrto = trch5pg.attrs.get('yrto')
yrcnt = yrto - yrfrom + 1
orglist = trch5pg.attrs.get('organs')
orgcnt = len(orglist)

# write attributes to new file
trch5fx.attrs.modify('nreps',nreps)
trch5fx.attrs.modify('npeople',npeople)
trch5fx.attrs.modify('yrfrom',yrfrom)
trch5fx.attrs.modify('yrto',yrto)
trch5fx.attrs.modify('organs',orglist)

trch5fx.attrs.items()

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

zdose = np.zeros((nreps,yrcnt,orgcnt))    
noxcnt = 0
noicnt = 0
inegcnt = 0
ineglist = []
xnegcnt = 0
xneglist = []
noxlist = []
noilist= []
pcnt = 0
print

for people in trcmembers:
  pcnt +=1
  if pcnt % 100 == 0:
      print pcnt,
      if pcnt % 1000 == 0:
          print

  try:
      pgrp = trch5pg[people]
      rzdids = pgrp['int'][:]
      # identify negative doses and set to 0      
      if np.min(rzdids) < 0:
          inegcnt +=1
          ineglist.append(people)
          rzdids[rzdids < 0] = 0
  except:
      print '\n No int', people, pcnt
      noicnt += 1
      noilist.append(people)
      dname = people + '/int'
      rzdids = zdose
      
  trch5fx.create_dataset(people+'/int', data=rzdids, compression='gzip', compression_opts=9)      
  # get external doses if they exist otherwise define them as 0
  try:
    rzdxds = pgrp['ext']
      # identify negative doses and set to 0      
    if np.min(rzdxds) < 0:
        xnegcnt +=1
        xneglist.append(people)  
        rzdxds[rzdxds < 0] = 0
  except:
    noxcnt += 1
    noxlist.append(people)
    dname = people + '/ext'
    rzdxds = zdose

  trch5fx.create_dataset(people+'/ext',data=rzdxds, compression='gzip', compression_opts=9)
     
print pcnt 

if noxcnt>0:
    if noxcnt == 1: 
        print noxcnt, "person with no external dose"
    else:
        print noxcnt, "people with no external dose"
    print noxlist

if xnegcnt>0:
    if xnegcnt == 1: 
        print '\n', xnegcnt, "person with negative external cum dose"
    else:
        print '\n', xnegcnt, "people with negative external cum dose"
    print xneglist

if noicnt>0:
    if noicnt == 1: 
        print noicnt, "person with no internal dose"
    else:
        print noicnt, "people with no internal dose"
    print noilist

if inegcnt>0:
    if inegcnt == 1: 
        print '\n', inegcnt, "person with negative external cum dose"
    else:
        print '\n', inegcnt, "people with negative external cum dose"
    print ineglist

trch5pg.close()
trch5fx.close()

 
        


