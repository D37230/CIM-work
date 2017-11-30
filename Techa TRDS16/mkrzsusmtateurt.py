# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 16:42:51 2017

@author: Dale
"""
import h5py
import os
import sys
import numpy as np


ghdir = 'C:/Users/prest/Documents/GitHub/CIM-Work/'
cimworkdir = ghdir + 'Techa TRDS16'
sys.path.append(cimworkdir)

import hdf5utils as trch5


trch5dir = 'h:/uraldosimetry/techa/hdf5/'

h5base = 'eurt16mc'
h5rz  = 'rz1'

rzname = trch5dir +h5base +h5rz+'.h5'
os.chdir(trch5dir)


trch5rz = h5py.File(rzname,'r')

idnumlist = trch5rz['sumstat/idlist'][:]





rzh5 = trch5rz['realizations']
ssh5 = trch5rz['sumstat']

[rzam, rzasd,rzgm, rzgsd] = trch5.rzsumstat(rzh5,'sto')
np.savetxt('eistoam.txt',rzam[:,:,0],delimiter=',')
np.savetxt('exstoam.txt',rzam[:,:,1],delimiter=',')
np.savetxt('eistogm.txt',rzgm[:,:,0],delimiter=',')
np.savetxt('exstogm.txt',rzgm[:,:,1],delimiter=',')
np.savetxt('eurtids.txt',idnumlist,fmt='%d')
xxx
for organ in orglist:
    [rzam, rzasd,rzgm, rzgsd] = trch5.rzsumstat(rzh5,organ)
    try:
        ssorg = ssh5.create_group(organ)
    except:
        print 'sumstat', organ, 'group already exists'
        ssorg = ssh5[organ]
    for stats in ['am', 'asd', 'gm', 'gsd']:
        dsname =  'rz'+stats
        print dsname
        # create statistics dataset in sumstat/organ group        
        try:
            ssorg.create_dataset(stats,data=dsname, dtype = 'f', compression = 'gzip',
                                compression_opts=9)
        except:
            print 'sumstat', organ, stats, 'dataset already exists'
    
trch5rz.close()
trch5pg.close()           
            