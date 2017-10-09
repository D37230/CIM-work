# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 07:11:32 2017

@author: prest
"""
import numpy as np
import pandas as pd
import h5py

def h5update(h5pg,popnm, batch, batchnm, udoses):
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
            
            udoses
                names of doses (columns to be included in hdf5 files)
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
    # print h5pg.keys()
    lastkeycnt = len(h5pg.keys())
    print lastkeycnt
    
    ihdr = pd.read_csv(icsvfn, nrows=0, sep='\s+')
    ihdr = ihdr.rename(columns = {'ubid':'sysnum'})
    
    # get dose column names and sert dose names
    dscols =  ihdr.columns.tolist()[3:] 
    dsnames  =  [dn.replace(iext,'') for dn in dscols]
    for newnm in dsnames:
        oldnm = iext + newnm
        ihdr = ihdr.rename(columns = {oldnm:newnm})
    
    ihdr = ihdr.rename(columns = {'uli':'col'})
    ihdr = ihdr.rename(columns = {'lli':'rec'})
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

    updatepg(h5pg, icsvfn, 'int', ihdr,dsnames, udoses, nreps, yrcnt, pperck, ckmax )
    updatepg(h5pg, xcsvfn, 'ext', ihdr,dsnames, udoses, nreps, yrcnt, pperck, ckmax )
    
    h5pg.attrs['yrfrom']=yrfrom    
    h5pg.attrs['yrto']=yrto
    h5pg.attrs['organs'] = udoses
    h5pg.attrs['units'] = 'Gy'
    h5pg.attrs['dlayers'] = ['int', 'ext']
    
    keylist = h5pg.keys()
    newkeycnt = len(keylist)
    newkeys = newkeycnt  - lastkeycnt
    print '\n',newkeys, "members added to ", popnm, " HDF5 from batch", batchnm
    print 'Total members in', popnm, 'HDF5 file:', newkeycnt
    
    np.savetxt(popnm+'keylist.csv',keylist,fmt='%s')

def updatepg (h5file, dosecsv, dtype, dshdr, dscols, udoses, nreps, yrcnt, pperck, ckmax ):
    ckno = 0
    for chunk in pd.read_csv(dosecsv, header=None, skiprows=[0],
                      sep='\s+', names = dshdr, chunksize=nreps*yrcnt*pperck, 
                      iterator=True, na_values='********'): 
        chunk.fillna(value=0, inplace=True)
        people = pd.unique(chunk['sysnum'])
        print ckno, len(chunk), people
        chunk = chunk.sort_values(by=['sysnum', 'real', 'year'])
        for p in people:
            addtopg(h5file, p, chunk, dtype, dscols, udoses,nreps, yrcnt)
        ckno += 1
        if ckno >= ckmax:
            break

def addtopg(h5pg, p, dsdat, dstype, dscols, udoses, rlzcnt, yrcnt, pu='', px=''):
    '''
    dosestoh5pg - save doses in csv file to HDF5 file with person groups
    
    Arguments:
            h5pg
                hdf5file for doses
            
            p
                person identifier

            dsdat
                pandas data frame with dose realization data
                (one row per person-year) with columns sysnum, real, year, and doses
            
            dstype
                'int' for internal or 'ext' for external
            
            dscols
                list of dose column names in input file
                
            udoses
                list of columns to be used
                

    '''
    pstr = str(p)
    doses = dsdat.loc[dsdat['sysnum']==p, 'real':dscols[-1]]
    isexposed = np.max(doses[dscols[1]])>0
    if isexposed:
        # fill 0 dose with previous years value
        doses = doses.apply(lambda x: x.replace(to_replace=0, method='ffill'), axis=0)
        try:
            dsarr = doses.as_matrix(columns=udoses)
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
            
    return 

        
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


def  pgtorlz(pgh5, rzh5, npeople, nreps, yrcnt, orgno, organ, rzrows, chunksz=50):
    '''
    pgtorlz  move for an organ data in a person-group hdf5 file to a relaization group hdf5 file
    
    Arguments:
        pgh5
            Person group hdf5 file
        
        rzh5
            Realization group hdf5 file (usually the 'realizations' group in
            an hdf5 file with 'realzations' and 'sumstat' groups)
            
        npeople
            Number of people in with dose
        nreps
            Number of realizations
            
        yrcnt
            NUmber of years
            
        orgno
            layer number for organ of interest
            
        organ
            name for organ of interest
            
        rzrows
            dictionary of Id's with link to row for that person in 
            realization datasets
            
        chunksz
            Number of relaizations to process per pass through people in 
            person group file
    
    Notes:
        
        person-group hdf5 structure:
            /ID groups with 3D int and ext datasets
            realization by year by organ
        
        realization-group structure
            /realizaitons/rzno with 3D organ datasets
            person by year by [int, ext]
            
    ''' 
    
    print '\n',organ,
    chunksz = 50
    nchunks = nreps/chunksz
    # initialize realization data arrayfor this organ
    rlzdata = np.zeros((chunksz,npeople,yrcnt,2))

    for chunk in range(nchunks):
        rznos = range(chunk*chunksz,(chunk+1)*chunksz)
        for rzno in rznos:
            rzstr = str(rzno)
            try:
                rzngrp = rzh5.create_group(rzstr)
            except:
                rzngrp = rzh5[rzstr]
            
        i = 0
        print '\nChunk ',chunk,':',
        for people in pgh5.iteritems():
            i += 1
            if i % 5000 == 0:
                print i,
            rlzrow = rzrows[people[0]]
            # print people, rlzrow, rzno, orgno
            rzdi = people[1]['int'][rznos,:,orgno] 
            rlzdata[range(0,chunksz),rlzrow,:,0] = rzdi
            rzdx = people[1]['ext'][rznos,:,orgno]
            rlzdata[range(0,chunksz),rlzrow,:,1] = rzdx

        print i, "Creating HDF5 datasets",
        for i in range(0,chunksz):
            if i % 10 == 0:
                print rznos[i],
            rzngrp = rzh5[str(rznos[i])]
            rzngrp.create_dataset(organ,data=rlzdata[i,:,:,:], dtype='f',
                                     compression='gzip', compression_opts=9)

        print rznos[i],    
# end of pgtorlz        