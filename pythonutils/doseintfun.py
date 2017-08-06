# -*- coding: utf-8 -*-
"""
Dopse interpolation functions:
    setdoseintinfo      Setupo for dose tnterolation
    doseinterp          Dose interp0olation (based on info from setdoseintindo)
Created on Thu Feb 11 10:31:44 2016

@author: Dale Prestonw
"""
import numpy as np


def setdoseintinfo(rowindx, minfo, idvar, yrvar, lag, fromyr, toyr):
    '''
    Setup for dose4 interpolation
    Written by Dale Preston
    February 2016

    Arguments:
        rowindx     associative array linking people to rows in person-year
                    dose matrix from HDF5 realization-hierarchy file
        minfo       data set with informaiton about ungrouped peson year
                    file -- must include named columns for id (numeric)
                    and fractional years
        idvar       name of ID variable column in minfo
        yrvar       name of year variable column in minfo
        lag         lag for dose computaiton (years)
        fromyr      first year in dose matrix
        toyr        last year in dose matrix

    Returns
        rowno       vector with dpse matrixrow number for each row in
                    minfo
        locol       dose-matrix column number for dose at start of
                    exposure year for each person
        hicol       dose-matrix column number for dose at end of
                    exposure year for each person
        uyr         fractional years from minfo
        ufrac       fraciton of year from uyr

    '''
    # define values for interpolation
    nyrs = toyr - fromyr + 2
    ids = minfo[idvar]
    idstr = [str(int(x)) for x in ids]
    rowno = np.zeros(len(idstr))
    for mrow, id in enumerate(idstr):
        try:
            rowno[mrow] = rowindx[id]
        except:
            rowno[mrow] = -1
    # rowno = [rowindx[str(x)] for x in minfo[idvar]]  # dose row number
    # adjusted fractional year for lookup
    uyr = [min(max(x-lag, fromyr), toyr) for x in minfo[yrvar]]
    # first year column
    locol = [max((int(x)-fromyr), 0) for x in uyr]
    # last year column
    hicol = [min(x+1, nyrs) for x in locol]
    # year fractions
    ufrac = [x - int(x) for x in uyr]

    intinfo = {}
    intinfo.update({'rowno' : rowno,
                    'loidx'  : locol,
                    'hiidx' : hicol,
                    'uyr' : uyr,
                    'ufrac' : ufrac})

    return(intinfo)


def doseinterp(doserlz, rlnostr, dorg, imdinfo, 
               scale=1, dlayer=-1):
    '''
    Dose interpolation function
    Written by Dale Preston
    Feburary 2016

    Extracts a person-by-year dose matrix for a specific
    "realizaiton" (rlnostr) of dose (dorg) from an HDF5 dose
    realization file (doserlz)
    Uses this data to interpolate doses for in a
    person year vector with fractional years.  The data needed to do this
    is summarized by 4 vectors:
    1)  a list that idetnifies the dose-matrix row (person) for each entry
       in the person-year matrix
    2) vector of dolse-matrix column numbers for the annnual dose at the start
       of the year of interest (locol)
    3) vector of dolse-matrix column numbers for the annnual dose at the start
       of the year of interest (hicol)
    4) vector of frational years (range 0 to 1) for the interpolation (ufrac)

    Arguments:
        doserlz  -  dose realziation hdf5 file with person-by=-year matrix for
                    each realizaiton
        rlnostr  -  realizaiton string.  This can be a string representation
                    of a realization number or `sumstat' if one wants mean
                    doses
        dorg     -  The HDF% group name for the organ of interest.  e.g. 'rbm'
                    or 'lng'.  For mean doses this should be a string with the
                    organ group name and 'am' as a subgroup, e.g. 'rbm/am'
        imdinfo   - interpolation data info (rownumbers, hi and low column indicies, fractional years)
        dlayer    - dose layer to use when dose matrix contains multiple
                    layers (e.g. external and internal doses.  This is
                    ignored if the dose arry has only two imensions.
                    The dfault value (-1) indicates that the doses should
                    be summed over the layers.
    '''
    rlzdnam = rlnostr + '/' + dorg
    # print rlzdnam
    # print doserlz
    usdose = doserlz[rlzdnam][:]
    # compute total doses if third dose dimension
    if usdose.ndim == 3:
        if dlayer < 0:
            usdose = np.sum(usdose, axis=2)
        else:
            usdose = usdose[:, :, dlayer]
    doses = [x*scale for x in usdose]

    # add initial column of zeors
    doses = np.c_[np.zeros(len(doses)), doses]
    # print doses.shape
    rowno = imdinfo['rowno']
    idose = np.zeros(len(rowno))
    for i in enumerate(rowno):
        drow = int(i[1])
        mrow = int(i[0])
        if drow >= 0:
            idose[mrow] = imdinfo['ufrac'][mrow] * (doses[drow, imdinfo['hiidx'][mrow]] -
                                         doses[drow, imdinfo['loidx'][mrow]])
            idose[mrow] += doses[drow, imdinfo['loidx'][mrow]]
        else:
            idose[mrow] = 0
    return np.reshape(idose,(-1,1))


def getcumtot(doserlz, relnostr, dorg, scale= 1):
    '''
    Get vector of total cumulative doses for a specified "realization"

    Arguments:
        doserlz     hdf5 realization file with dose person-year dose arrays
        rlznostr    realization string  can be either a realization number or
                    'susmtat' for the means
        dorg        organ string 'org/dose' for realizaiton info or
                    'org/am'  or 'org/gm'  for sumstat info

    Returns
       cumdose     vector of cumdoses
    '''
    cd = np.max(doserlz[relnostr + '/' + dorg][:], axis=1)
    cumdose = [x*scale for x in cd]
    return cumdose


def showintinfo(doserlz, pindx, dnam, id, fromyr, ugdose, rowno,
                uyr, ufrac, locol, hicol):
    idstr = str(id)
    xd = np.column_stack((np.zeros(len(pindx)), doserlz[dnam][:]))
    dhist = doserlz[dnam][pindx[idstr], :]
    print '\nDose history ID ', id, '\n'
    for i, x in enumerate(dhist):
        print fromyr+i, x
    print '\n Interpoalted dose info\n'
    condl = [x == pindx[idstr] for x in rowno]
    for i, x in enumerate(condl):
        if x:
            xlo = xd[pindx[idstr], locol[i]]
            xhi = xd[pindx[idstr], hicol[i]]
            ixd = xlo + ufrac[i] * (xhi - xlo)
            print uyr[i], ufrac[i], locol[i], hicol[i], xlo, xhi, \
                ixd, ugdose[i]
    print '\n'




    