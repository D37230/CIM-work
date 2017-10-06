# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 15:25:55 2017

@author: Dale
"""

import h5py
import os
import numpy as np

'''
  Fix zedro dose problems in TRC person group files
  a)  initial zero doses for realizations after the first one were incorrectly filled
      with value of total cumulative dose for previous realization
  b)  datasets for 260 TRC members with non-zero internal doses but 0 external doses 
      were not written to the hdf5 file
  c)  individual with 0 internal and external doses was not included in file (sysnum 132119)

'''

trch5dir = 'k:/uralsdosimetry/techa/hdf5'
h5base = 'trc16mc'
h5sub = 'ds1'

os.chdir(trch5dir)

trch5pg = h5py.File(h5base+h5sub+'.h5','r')
trch5fx = h5py.File("x:/"+h5base+h5sub+'fix.h5')


    
trcmembers = trch5pg.keys()

nreps = 1500
npeople = len(trcmembers)
yrfrom = trch5pg.attrs.get('yrfrom')
yrto = trch5pg.attrs.get('yrto')
yrcnt = yrto - yrfrom + 1
orglist = trch5pg.attrs.get('organs')
orgcnt = len(orglist)


# create sorted list of id's for realization data set
idlist = []
for chid in trcmembers:
    idlist.append(int(chid))

zdose = np.zeros((nreps,yrcnt,orgcnt))    
noxcnt = 0
negcnt = 0
neglist = []
noxlist = []
pcnt = 0
print
for people in trch5pg:
  pcnt +=1
  if pcnt % 100 == 0:
      print pcnt,
  pgrpfx = trch5fx.create_group(people)
  rzdi = trch5pg[people + '/int'][:]
  # get external doses if they exist otherwise define them as 0
  try:
    rzdx = trch5pg[people + '/ext'][:]
  except:
    noxcnt += 1
    noxlist.append(people)
    rzdx = zdose
  if (np.min(rzdi) < 0 ) | (np.min(rzdx) < 0 ):
      negcnt += 1
      neglist.append(people)

  # check the number of leading nonzero doses for the first realization
  # set the leading elements to 0 for all realizations
  r0 = rzdi[0,:,0]
  lzcnt = 0
  for yr in range(0,yrcnt):
    if(r0[yr] == 0):
      lzcnt += 1
    else:
      break
  if lzcnt > 0:
    rzdi[:,range(lzcnt),:] = 0
    rzdx[:,range(lzcnt),:] = 0

  # put corrected internal and external datasets in to the new hdf5 file 
  pgrpfx.create_dataset("int", data=rzdi, compression='gzip', compression_opts=9)
  pgrpfx.create_dataset("ext", data=rzdx, compression='gzip', compression_opts=9)

# add group for person with 0 internal and external MC dose esitmates
pgrpfx = trch5fx.create_group('132119')
pgrpfx.create_dataset("int", data=zdose, compression='gzip', compression_opts=9)
pgrpfx.create_dataset("ext", data=zdose, compression='gzip', compression_opts=9)


if noxcnt>0:
    if noxcnt == 1: 
        print noxcnt, "person with no external dose"
    else:
        print noxcnt, "people with no external dose"
    print noxlist

if negcnt>0:
    if negcnt == 1: 
        print '\n', negcnt, "person with negative cum dose"
    else:
        print '\n', negcnt, "people with negative cum dose"
    print neglist


trch5fx.close()


 
        


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