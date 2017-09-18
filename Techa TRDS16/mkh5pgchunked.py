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

sys.path.append('c:/Users/prest/Documents/GitHub/CIM-work/Techa TRDS16')

import hdf5utils as trch5


dsrlzindir = "h:/Uraldosimetry/techa/napier/"
dsrlzindir = "e:/"
dsrlzh5dir = "e:/hdf5/"

# from hdf5utils import addtoh5pg


os.chdir(dsrlzh5dir)



# batchstats = trch5.batchinfo(icsvfn, iext)

def trch5update(h5fn,popnm, batch, batchnm):
    '''
        trch5update :  Update TRC/EURT hdf5 file using raw TRDS16 output files
        
        Arguments:
            
            h5fn
                HDF5 file name  (opened and closed inside this function)
            
            popnm
                population name:  'trc' or 'eurt'
            
            batch
                batch file suffix (like '2k-3k')
                
            batchnm
                root of batchfile name with path but no extension
                extension (ein/tin, eou/tou) added in this routine
    '''
    if popnm == "trc":
        iext = "tin"
        xext = "tou"
    elif popnm == "eurt":
        iext = "ein"
        xext = "eou" 
        
    icsvfn  = batchnm + "." + iext
    xcsvfn  = batchnm + "." + xext
    print icsvfn
    h5pg = h5py.File(h5fn)
    lastkeycnt = len(h5pg.keys())
    print lastkeycnt
    
    ihdr = pd.read_csv(icsvfn, nrows=0, sep='\s+')
    ihdr = ihdr.rename(columns = {'ubid':'sysnum'})
    
    # get dose column names and sert dose names
    dscols =  ihdr.columns.tolist()[3:] 
    dsnames  =  [dn.replace(iext,'') for dn in dscols] 
    # dscnt = len(dsnames)
    
    nreps = 1500 
    yrfrom = 1950
    if popnm == "eurt":
        yrfrom = 1957
    yrto = 2015
    yrcnt = yrto-yrfrom+1
    
    # set up number of chunks
    pperck = 20
    ckmax = 1000/pperck
    if ckmax*pperck < 1000:
        ckmax +=1
       
    ckno=0
    
    # h5pgupdate (h5file, dosecsv, dshdr, dscols, nreps, yrcnt, pperck, ckmax ):
    trch5.updatepg(h5pg, icsvfn, 'int', ihdr,dscols, nreps, yrcnt, pperck, ckmax )
    trch5.updatepg(h5pg, xcsvfn, 'ext', ihdr,dscols, nreps, yrcnt, pperck, ckmax )
    
    
    h5pg.attrs['yrfrom']=yrfrom    
    h5pg.attrs['yrto']=yrto
    h5pg.attrs['organs'] = dsnames
    h5pg.attrs['units'] = 'Gy'
    h5pg.attrs['dlayers'] = ['int', 'ext']
    
    keylist = h5pg.keys()
    newkeycnt = len(keylist)
    newkeys = newkeycnt  - lastkeycnt
    print '\n',newkeys, "members added to ", popnm, " HDF5 from batch", batchid
    print 'Total members in', popnm, 'HDF5 file:', newkeycnt
    
    np.savetxt(popnm+'keylist.csv',keylist,fmt='%s')
       
    h5pg.close()


# TRC HP5 person group file
batchpfx = "2017Cases"
batchids = ['1-1000', '1k-2k', '2k-3k', '3k-4k','4k-5k', '5k-6k', '6k-7k', '7k-8k',
            '8k-9k', '9K-10K',  '10k-11k', '11K-12K', '12k-13k', '13k-14k' '14k-15k',
            '15k-16k', '16k-17k', '17k-18k']


cohort = "trc"
h5filenm = dsrlzh5dir + cohort + 'x16mc.h5'
batchlist = ['12k-13k', '13k-14k']
for batchid in batchlist:
    batchnm  = dsrlzindir + batchpfx + batchid 
    print batchid
    trch5update(h5filenm,cohort,batchid,batchnm) 


cohort = "eurt"
h5filenm = dsrlzh5dir + '/x' + cohort + '16mc.h5'
['15k-16k', '17k-18k']
for batchid in batchlist:
    batchnm  = dsrlzindir + batchpfx + batchid 
    print batchid
    trch5update(h5filenm,cohort,batchid,batchnm)


