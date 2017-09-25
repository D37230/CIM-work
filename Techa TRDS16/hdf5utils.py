# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 07:11:32 2017

@author: prest
"""
import numpy as np
import pandas as pd
import h5py

def addtopg(h5pg, p, dsdat, dstype, dscols, rlzcnt, yrcnt, pu='', px=''):
    '''
    dosestoh5pg - save doses in csv file to HDF5 file with person groups
    
    Arguments:
            h5pg
                hdf5file for doses
            
            p
                person identifier

            dsdat
                pandas data fream with dose realization data
                (one row per person-year) with columns sysnum, real, year, and doses
            
            dstype
                'int' for internal or 'ext' for external
            
            dscols
                list of dose column names
                
            pu
                string of updated data tables (int, ext dsmn)
                
            px
                string of no updated data tables
    '''
    pstr = str(p)
    doses = dsdat.loc[dsdat['sysnum']==p, 'real':dscols[-1]]
    dsmn=0
    isexposed = np.max(doses[dscols[0]])>0
    if isexposed:
        # fill 0 dose with previous years value
        doses = doses.apply(lambda x: x.replace(to_replace=0, method='ffill'), axis=0)
        dsmn = 0 # np.mean(dsarr,axis=0)
        try:
            dsarr = doses.as_matrix(columns=dscols)
            dsarr = dsarr.reshape((rlzcnt,yrcnt,-1),order='C')
            pgroup = pstr + '/' + dstype
            try:
                h5pg[pgroup][:]
                px += dstype + ' '
            except:
                pu += dstype + ' ' 
                h5pg.create_dataset(pgroup,data=dsarr,dtype='f',compression='gzip',compression_opts=9)
        except:
            print "***** Bad data ***** data sysnum ",p
            
    return pu, px, dsmn

def updatepg (h5file, dosecsv, dtype, dshdr, dscols, nreps, yrcnt, pperck, ckmax ):
    ckno = 0
    for chunk in pd.read_csv(dosecsv, header=None, skiprows=[0], 
                      sep='\s+', names = dshdr, chunksize=nreps*yrcnt*pperck, 
                      iterator=True, na_values='********'): 
        chunk.fillna(value=0, inplace=True)
        people = pd.unique(chunk['sysnum'])
        print ckno, len(chunk), people
        chunk = chunk.sort_values(by=['sysnum', 'real', 'year'])
        for p in people:
            pu, px, idsmn = addtopg(h5file, p, chunk, dtype, dscols, 
                                      nreps, yrcnt)
        ckno += 1
        if ckno >= ckmax:
            break
        
def batchinfo(batchfn, iext, chksz = 300000, usecols= 3):
    '''
    bathcinfo: get basic info about TRDS16 MC batch 
    
    Arguments:
        batchfn
            file name for batch (like 2017cases1k-2k)
        iext
            file name extension (assumed prefix for organ doses in column names)
        chksz
            number of rows to read_csvd for checkinf number of realizations and 
            year range.  Default 300,000
        usecols
            number of columns to read. Default 3 (sysnum, relizaiton number, year)

    Returns:
        dictionary with entries
        
            dscols: dose column cname in pandas dataframe
            
            dsname: organ names
            
            yrfrom: first year
            
            yrto:   last year
            
            nreps: number of dose realizatons
       
    '''
    hdr = pd.read_csv(batchfn, nrows=0, sep='\s+')
    hdr = hdr.rename(columns = {'ubid':'sysnum'})
    dscols =  hdr.columns.tolist()[3:] 
    dsnames  =  [dn.replace(iext,'') for dn in dscols] 
    dscnt = len(dsnames)
    finfo = pd.read_csv(batchfn, header=None, skiprows=[0],usecols=range(usecols),
                             sep='\s+', names = hdr, nrows = chksz)
    idcnt = len(pd.unique(finfo['sysnum']))
    yrcnt = len(pd.unique(finfo['year']))
    yrfrom = np.min(finfo['year'])
    yrto = np.max(finfo['year'])
    rlzcnt = len(pd.unique(finfo['real']))
    print '\n', rlzcnt, "realizations with",yrcnt, 'years from ', yrfrom, 
    print 'to',yrto,'for', dscnt, 'organ doses' ,'\n'
    print dsnames
    runinfo = {}
    runinfo.update({'dscols': dscols,
                    'dsnames': dsnames,
                    'yrfrom': yrfrom,
                    'yrto' : yrto,
                    'yrcnt': yrcnt,
                    'nreps': rlzcnt})
    return runinfo

def h5update(h5fn,popnm, batch, batchnm):
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
    h5pg = h5py.File(h5fn,'a')
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
    
    # h5pgupdate (h5file, dosecsv, dshdr, dscols, nreps, yrcnt, pperck, ckmax ):
    updatepg(h5pg, icsvfn, 'int', ihdr,dscols, nreps, yrcnt, pperck, ckmax )
    updatepg(h5pg, xcsvfn, 'ext', ihdr,dscols, nreps, yrcnt, pperck, ckmax )
    
    
    h5pg.attrs['yrfrom']=yrfrom    
    h5pg.attrs['yrto']=yrto
    h5pg.attrs['organs'] = dsnames
    h5pg.attrs['units'] = 'Gy'
    h5pg.attrs['dlayers'] = ['int', 'ext']
    
    keylist = h5pg.keys()
    newkeycnt = len(keylist)
    newkeys = newkeycnt  - lastkeycnt
    print '\n',newkeys, "members added to ", popnm, " HDF5 from batch", batchnm
    print 'Total members in', popnm, 'HDF5 file:', newkeycnt
    
    np.savetxt(popnm+'keylist.csv',keylist,fmt='%s')
       
    h5pg.close()    
        