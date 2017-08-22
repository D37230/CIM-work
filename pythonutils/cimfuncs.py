# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 12:41:37 2016

@author: Dale Preston
"""
import numpy as np
from scipy.stats import lognorm, norm
from scipy.integrate  import quad
from scipy.optimize  import brentq
from doseintfun import doseinterp

def mkQmat(stcovs, fv, mfit, modtype):
    '''
    Define Q matrix    $d log(\mu) / d \beta
    
    Currently only for mulitplcative ERR Model 
    but easily modified tto compute for multiplicative ERR additive ERR and EAR
    based on modtype argument
    
    Arguments:
        
        stcovs 
                list of subeterm specific covariate matrices
        
        mfit 
                model fit structure based on psv info as returned from modeleval
                computations use term values ('termvals' and subterm values ('subtvals')\
                
        mu
            vector of fitted values
        
        modtype 
                model type from  psv file obtained form 'type' element of 
                the model fit struture  
                
                1 - Additive ERR (T_0 (1 + T_1 + T-2 + ...)
                
                3 - EAR T0 + T1 + T2 + ... 
                
                5 - Multiplicative ERR To ( 1 + T1) (1 + T2) ...
                  
    Returns:
        
        Q - p x p matrix 
        
    '''
    Q = ()
    multerr = 0 + (modtype == 5)
    adderr = 0 + (modtype == 1)
    earmod = 0 + (modtype == 3)
    back = np.ones_like(mfit['termvals'][0][1])
    t0 = mfit['termvals'][0][1]
    if adderr : 
        back = t0  # assumes that there is a T0
    elif earmod:
        back = np.ones_like(t0)
    for stno in range(len(stcovs)):
        tno = mfit['subtvals'][stno][0]     # termno
        isloglin = mfit['subtvals'][stno][1]==1
        stval = mfit['subtvals'][stno][2]
        tval = mfit['termvals'][tno][1]
        tfit = tval
        if(tno>0):
            tfit = t0*tval       
 
        trem = np.ones_like(tval)
        trem[stval!=0] = tval[stval != 0] /stval[stval !=0]
        # termerr = mu/t0
        if multerr:
            back = fv/((tno>0) + tval) #  mu/((tno>0) + tval)
        
        # temp = stcovs[stno]/(cons  + mfit['termvals'][tno][1])
        temp = stcovs[stno] * trem * back
        # print temp.shape
        if isloglin :    # logl linear subterm
            temp = temp * stval
        if stno == 2:
            # temp = stcovs[2] * tval/t0 * mu / (1 + tval/t0)
            temp = temp
            #  print "term ", mfit['subtvals'][stno][0], isloglin, np.shape(temp)
        if stno == 0:
            Q =  temp
        else:
            Q = np.column_stack((Q,temp))
    return Q/fv


def mkcim(Q, G, rlzdb, rlzcnt, dose, interpinfo, 
            prmcov, dslayer=0,
          infoevery=100, scale=1):
    '''  
    Compute CIM adjustment matrix

    Arguments:

        Q
            Q matrix (score like as per computation notes)

        G
            Diagonal of G matrix (as per computation notes)

        rlzdb 
            HDF5 dose realization file
            
        rlzcnt
            number of realizations

        interpinfo
            basic info for interpolation (created using setdoseintinfo)
        
    Returns:
        
        CIM matrix
            I-1 M I-1
    '''
    QG = np.transpose(Q*G)
    rlzno = 0
    infotimes = 0
    print "Number of realizations ", rlzcnt
    print QG.shape
    for x in range(rlzcnt):
        idrlz = doseinterp(rlzdb, str(int(x)), dose, interpinfo,
                           dlayer=dslayer, scale=scale)
        idrlz = np.reshape(idrlz, (-1, 1))
        if x == 0:
            qgd = np.dot(QG, (idrlz))
        else:
            qgd = np.column_stack((qgd, np.dot(QG, idrlz)))
        rlzno += 1
        if infoevery:
            if (rlzno % infoevery) == 0:
                infotimes += 1
                print rlzno, " ",
                if infotimes == 10:
                    infotimes = 0
    cim = np.cov(qgd)
    return np.dot(prmcov, np.dot(cim, prmcov))

def mkcimdr(Q, G, intdreps, rlzcnt, prmcov, reppfx= 'rep', infoevery=100):
    '''  
    Compute CIM adjustment matrix

    Arguments:

        Q
            Q matrix (score like as per computation notes)

        G
            Diagonal of G matrix (as per computation notes)

           
        rlzcnt
            number of realizations

        prmcov
            parameter covariance matrix
            
        reppfx
            Replicaiton variable name prefix (defualt 'rep') names are
            assumed to have the form reppfxn  with n is the replication number
            
        infoevery
            print porgress indicator every infoevery realizations
        
    Returns:
        
        CIM matrix
            I-1 M I-1
    '''
    drows = intdreps['rowno'].astype(int)
    drep = np.zeros_like(G)
    QG = np.transpose(Q*G)
    rlzno = 0
    infotimes = 0
    print "Number of realizations ", rlzcnt
    for rno in range(rlzcnt):
        repname = reppfx + str(rno)
        idrlz = intdreps[repname]
        idrlz = np.reshape(idrlz, (-1, 1))
        drep[drows] = idrlz
        if rno == 0:
            qgd = np.dot(QG, (drep))
        else:
            qgd = np.column_stack((qgd, np.dot(QG, drep)))
        rlzno += 1
        if infoevery:
            if (rlzno % infoevery) == 0:
                infotimes += 1
                print rlzno, " ",
                if infotimes == 10:
                    infotimes = 0
    cim = np.cov(qgd)
    return np.dot(prmcov, np.dot(cim, prmcov))

def logNormal2Normal(amean, asd):
    '''
    lnmusig - compute log-scale mean and sd given aritimetic mean and sd

    Arguments:

    amean -  arithmetic mean

    asd   -  arithmetic standard deviation
    '''
    lnvar = np.log(1 + pow(asd/amean, 2))
    lnmean = np.log(amean) - lnvar/2

    return (lnmean, np.sqrt(lnvar))


def plnorm(x, lnmn, lnse, lower_tail=1):
    if lower_tail:
        return lognorm.cdf(x, lnse, scale=np.exp(lnmn))
    else:
        return 1-lognorm.cdf(x, lnse, scale=np.exp(lnmn))


def qlnorm(x, lnmn, lnse):
        return lognorm.ppf(x, lnse, scale=np.exp(lnmn))


def dnorm(x, mn, se):
    return norm.pdf(x, mn, se)


def qnorm(x, mn, se):
    return norm.ppf(x, mn, se)


def getbndrange(bhat, bse, mse, alpha):
    """
    Purpose: get hi and low n=bound search limits

    Arguments:

    bhat  - parameter estimate

    bxse  - parameter standard error

    mse   - dose error parameter variance from M'VM

    alpha  - confidence level e.g. 0.5 for 95% CI

    Return
    b_lo, b_hi  - search range limits
    lmu, lse    - mean and variance of lognormal(0,mse)

    """
    lmu, lsd = logNormal2Normal(1, mse)
    if (bhat - qnorm(1-alpha/2, 0, bse)) < 0:
        b_lo = (bhat - qnorm(1-alpha/3, 0, bse)) / qlnorm(alpha/3, lmu, lsd)
    else:
        b_lo = (bhat - qnorm(1-alpha/3, 0, bse)) / qlnorm(1-alpha/3, lmu, lsd)
    if (bhat - qnorm(alpha/2, 0, bse) < 0):
        b_hi = (bhat - qnorm(alpha/3, 0, bse)) / qlnorm(1-alpha/3, lmu, lsd)
    else:
        b_hi = (bhat - qnorm(alpha/3, 0, bse)) / qlnorm(alpha/3, lmu, lsd)

    return (b_lo, b_hi)


def intfunpos(y, bhat, bse, b, lnmn, lnse):
    # print 'pos ', bhat, y, b, lnmn, lnse
    return plnorm((bhat-y)/b, lnmn, lnse) * dnorm(y, 0, bse)


def intfunneg(y, bhat, bse, b, lnmn, lnse):
    #  print 'neg ', bhat, y, b, lnmn, lnse
    return plnorm((bhat-y)/b, lnmn, lnse, 0) * dnorm(y, 0, bse)


def cdf_approx(b, bhat, bse, mse, offset=0):
    """
    Evaluate lower tail probability for normal-lognormal mixture
    for dose uncertainty adjustment

    Arguments:

    b     - parameter value for which cdf is desired

    bhat  - parameter estimate

    bse   - parameter estimate stbadard error

    mse   - MVM standard error for parameter Wald bounds returned if mse=0)

    offset - value to subtract from cdf value (default 0)

    Returns:

    (blower, bupper)   -  upper and lower bound estimates

    """
    lmu, lsd = logNormal2Normal(1, mse)
    if bse < 1e-6:
        if b > 0:
            return plnorm(bhat/b, lmu, lsd) - offset
        else:
            return plnorm(bhat/b, lmu, lsd, 0) - offset
    else:
        if b < 0:
            return quad(intfunpos, -np.inf, np.inf,
                        args=(bhat, bse, b, lmu, lsd))[0] - offset
        else:
            return quad(intfunneg, -np.inf, np.inf,
                        args=(bhat, bse, b, lmu, lsd))[0] - offset


def getCI(bhat, bse, mse, alpha):
    if mse == 0:
        blower = qnorm(alpha/2, bhat, bse)
        bupper = qnorm(1-alpha/2, bhat, bse)
    else:
        b_lo, b_hi = getbndrange(bhat, bse, mse, alpha)
        blower = brentq(cdf_approx, b_lo, b_hi,
                        args=(bhat, bse, mse, alpha/2))
        bupper = brentq(cdf_approx, b_lo, b_hi,
                        args=(bhat, bse, mse, 1-alpha/2))
    return (blower, bupper)

def dispcimbnds (desc, prmest, waldbnds, cimbnds, estfmt = '%6.3f'):
    '''
    Print wald abd CIM bounds for a parameter
    
    Arguments:
        
        desc
            Paramter description / name
            
        prmest
            Parameter estimate
            
        waldbnds
            Wald bound array
            
        cimbnds
            CIM-adjusted boudn array
            
        estfmt
            Output format for estimate and bounds
            
    Returns:
        None
    '''
    strfmt = '%' + str(max(len(desc),16)) + 's'
    print strfmt % desc, " ", estfmt % prmest
    print strfmt % "Wald bounds:","(", estfmt % waldbnds[0], ',', estfmt % waldbnds[1] , ')'
    print strfmt %  "Adjusted bounds:","(", estfmt % cimbnds[0], ',',  estfmt % cimbnds[1] , ')'


def showcimbnds(pno, prmvec, prmcov, cim, desc, cilev=0.05, expo=0):
    """
    Compute and diisplay CIM-adjusted and Wald bounds

    Arguments:
        pno      paameter number
        prmvec   parameter vector

        ormcov   parameter covariance matrix

        cim      Correct information matrix

        desc     parameter description

        cilev    confidence tail prob  (default 0.05 i.e. 95% CI)

        expo     non-zero to exponentiate estimate and bounds


    """
    pest = prmvec[pno]
    mse = np.sqrt(prmcov[pno, pno])
    pmse = np.sqrt(cim[pno, pno])

    pdadjbnd = getCI(pest, mse, pmse, cilev)
    pdwldbnd = getCI(pest, mse, 0, cilev)

    if(expo):
        pest = np.exp(pest)
        pdadjbnd = np.exp(pdadjbnd)
        pdwldbnd = np.exp(pdwldbnd)

    print "\n", desc, "estimate: ", pest
    print "Wald bounds:      ", pdwldbnd
    print "Adjusted bounds:  ", pdadjbnd

def modcimprmbnds(modelinfo, cim, level=95):
    '''
    Compute CIM adjsuted bounds for all parameters in the model
    
    Arguments:
    
        modelinfo
            Basic model informaiton structure obtained from psv file
            using getmodinfo (psvdict2)
            
        cim
            CIM adjustment matrix (from mkcim (cimfuncs))
            
        level
            Confidence level (default 95%)
    '''
    tprob = 1 - level/100.0
    nparms = len(modelinfo['prmvec'])
    prmse = np.sqrt(np.diag(modelinfo['prmcov']))
    isfree = modelinfo['prmstat'] == 0
    parms = modelinfo['prmvec']
    print "\n",level,"% bounds"
    for pno in range(nparms):
        parname = '\n' + str(pno+1) + ' ' + modelinfo['prmnames'][pno]        
        if isfree[pno]:
            mprmest = parms[pno]
            mprmse = prmse[pno]
            mmse = np.sqrt(cim[pno, pno])
            mpucim = getCI(mprmest, mprmse, mmse, tprob)
            mpuwald = getCI(mprmest, mprmse, 0,tprob)

            dispcimbnds(parname,mprmest,mpuwald, mpucim)
        else:
            print parname, "fixed/aliased"