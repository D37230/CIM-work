#!/usr/bin/env python
import re
import numpy as np

# first we need to identify where each model
# is located within the file
def getmodinfo(psvfile):
    """
    **Parse contesnts of epicure psv file saving contents in a list of dictionaries**\n 
    **input**:\n
    psvfile  - psv file handle\n
    **returns**:\n
    models   - list of dictionaries for each model in the file\n\n\n
    **Model dictionary keys**\n
    title    - user-defined name for model\n
    prog     - program in which model fit (AMFIT, GMBO, PEANUTS)\n
    type     - model type (form)\n
                1  Additive ERR\n
                2  Additive RR\n
                3  EAR\n
                4  Geometric mixture\n
                5  Mulitplicative ERR\n
    totprm   - total number of model parameters (not counting strata)\n
    freprm   - total number of free parameters \n
    strprm   - number of strata\n
    df       - degrees of freeom (AMFIT, unconditional GNBO) or number of risk sets (PEANUTS, conditional GMBO)\n
    dev      - deviance (AMFIT, unconditional GNBO) or -2*loglik (PEANUTS, conditional GMBO))\n
    prmvec   - parameter vector (totprm,1)\n
    prmcov   - parameter covariance matrix (totprm,totprm)\n
    prmstat  - parameter status (0 free; 1 fixed, not aliased; 2 aliased)\n
    prmnames - covariate names\n
    """
    
    indx = []
    # establish a regular expression
    models_pattern = "AMFIT|GMBO|PECAN|PEANUTS"
    models_regex = re.compile(models_pattern, re.IGNORECASE)

    # read PSV file contents into a lis buffer for processing 
    psave_content = psvfile.readlines()
    
    # determine initial lines for each model in psv file
    for i, line in enumerate(psave_content):
        if models_regex.search(line):
            indx.append(i)

    # loop through models creating basic info for model dictionaries
    models = []
    for x in indx:
        # tmp buffer
        model = {}
        # use the index (line number) to return the line
        # with the model information
        model_rec = psave_content[x].strip().split(',')
        # line index
        model = {'index': x}
        # model title (previous line)
        model.update({'title': psave_content[x-1].strip()})
        # program name
        model.update({'prog': model_rec[0]})
        # model code
        model.update({'type': int(model_rec[1])})

        # print model

        # total number of parameters
        model.update({'totprm': int(model_rec[2])})
        # number of free parameters
        model.update({'freprm': int(model_rec[3])})
        # number of records to read
        model.update({'userec': int(model_rec[4])})
        # number of model strata
        model.update({'strprm': int(model_rec[5])})
        # degrees of freedom
        model.update({'df': int(model_rec[6])})
        # deviance
        model.update({'dev': float(model_rec[7])})
        # add the model to the models buffer
        models.append(model)

    # loop throught models to get parameter and covaraite info
    for model in models:
        # set up our interval for the read
        # bgn == start line
        # end == end line
        bgn = model['index'] + 1
        end = model['totprm'] + bgn

        # these will be terms for each model
        subterminfo = []
        # parameter names
        pnams = []
        # free or fixed
        pfixd = []
        # parameter estimates for each model
        pestm = []
        # covariance matrix
        covmr = []
        
        # read the data
        trmno = -1
        subtno = -1
        stprmcnt = 0
        pno = -1
        stinfo = ()
        for i in range(bgn, end):
            pno += 1

            line = psave_content[i].strip().split(',')
            # replace epicure .
            line = [re.sub(r'^\.', '0', i) for i in line]
            # parameter name
            pnams.append(line[3].strip('"'))
            # fixed or free -- 0,free 1,fixed 2,redundant (fixed)
            pfixd.append(int(line[4]))
            # parameter estimate (to vector)
            # print line[5]
            pestm.append(float(line[5]))
            # covariance matrix fow
            covmr.append([float(i) for i in line[7:]])

            # create subterm info list
            # termno, subtermno, start parm, parameter count
            thistrm = int(line[1])
            thissubt = int(line[2])

            if ((thistrm != trmno) | (thissubt != subtno) | (i == end)):
                if(trmno >= 0):
                    stinfo.append(stprmcnt) # + (pno == (model['totprm']-1)))
                    stpv = pestm[stinfo[2]:stinfo[2]+stinfo[3]]
                    # print stpv
                    stinfo.append(stpv)
                    subterminfo.append(stinfo)

                stinfo = [thistrm, thissubt, pno]
                trmno = thistrm
                subtno = thissubt
                stprmcnt = 0

            stprmcnt += 1

        stinfo.append(stprmcnt)
        stpv = pestm[stinfo[2]:stinfo[2]+stinfo[3]]
        # print stpv
        stinfo.append(stpv)
        subterminfo.append(stinfo)
        # append/extend to the model buffer
        model.update({'subterms': subterminfo,
                      'prmnames': np.array(pnams),
                      'prmstat': np.array(pfixd),
                      'prmvec': np.array(pestm),
                      'prmcov': np.array(covmr)})

    return(models)
    # end of getmodinfo function


def showmodinfo(model, detail=0, showcov=0):
    """
      **Show model information from dictionary**
       Input:

     * model  - model info dictionary
     * detail - 0 basic info only
              1 show parameter vector info
     * showcov  0 no
              1 yes
    """
    print '\n', model['title']
    modtype = {1: 'Additive ERR', 2: 'Additive RR', 3: 'EAR', 4: 'Geometric mixture', 5: 'Multiplciative ERR' }
    print modtype[model['type']],' fit using ',model['prog']
    print model['totprm'], ' parameters (' , (model['freprm']-model['strprm']), ' free parameters)'
    print 'Deviance: ', model['dev']
    if(detail):
        print '\nParameter estimates'
        for i in range(0, model['totprm']):
            print (i+1), model['prmnames'][i], model['prmvec'][i], np.sqrt(model['prmcov'][i][i])
# end of shomodindo


def modeleval(modinfo, modcovs, pyoffset):
    """
    Evaluate model subterms, terms, and fitted values

    Arguments:
        modinfo     Model structure info from PSV file
        modcovs     List of covariate matrices defined from Ungrouped file'
        pyoffset    PYR value for computation of fitted values

    """
    # Consistency checks
    stinfo = modinfo['subterms']
    stcnt = len(stinfo)
    if stcnt != len(modcovs):
        print "Wrong number of subterms in covariate matrix list: "
        print "expected", stcnt, "but found ", len(modcovs)
        print 
        return(0)
    else:
        for i in range(0, stcnt):
            # print len(stinfo[i][4]), modcovs[i].shape[1]
            if len(stinfo[i][4]) != modcovs[i].shape[1]:
                print "Term ", stinfo[i][0], " Subterm ",
                print stinfo[i][1], "covariate matrix size error"
                print len(stinfo[i][4]), modcovs[i].shape[1]
                return(0)
    
    # print "Model OK"
    tno = -1
    termno = -1
    stvals = []
    trmvals = []
    trmval = np.ones_like(pyoffset)

    earmod = modinfo['type'] == 3
    for stno in range(0, stcnt):
        if tno != stinfo[stno][0]:
            tno = tno + 1
            tv = []
            if tno > 0:
                tv.append(termno)
                if(termno == 0 | earmod):
                    trmval = pyoffset*trmval

                tv.append(np.reshape(trmval,(-1,1)))
                trmvals.append(tv)
                trmval = 1
        stv = []
        termno = stinfo[stno][0]
        subtype = stinfo[stno][1]
        stv.append(termno)
        stv.append(subtype)
        stfv = np.dot(modcovs[stno], stinfo[stno][4])
        if subtype == 1:
            stfv = np.exp(stfv)
        stv.append(np.reshape(stfv,(-1,1)))
        stvals.append(stv)
        trmval = trmval * stfv

    tv = []
    tv.append(termno)
    tv.append(np.reshape(trmval,(-1,1)))
    trmvals.append(tv)

    terminfo = {}

    terminfo.update({'subtvals': stvals,
                     'termvals': trmvals})
    return(terminfo)
# end of modeleval


def calc_fvalues(modinfo, trmvals):
    """
        Compute fitted values: total, baseline, total excess, and
        subterm main effect excess
        
        Arguments:
        
        modinfo  -  PSV model structure info reated by psvdict.getmidinfo
        trmvals  -  term and subterm value list for this model created by 
                    psvdict.model_eval using the model info and a dataset
                    
        Return:
            List of fitted value structures with elements
            
            fv      -  total fitted value vector (nrec x 1)
            extot   -  total excess case vector  (nrec x 1)
                       (excess cases associated with  all non-baseline subterms)
            termfv  -  list of term-specific fitted case estimate vectors
                        one element for each term  elements [trmno,  trmval]
                        trmno   - termno 
                        trmval  - vector of term fitted cases (nrec x 1)
            exjnt   -  joint effect excess (non-zero only for multiplicative ERR only)
                       difference between total excess and sum of term-specific excess cases
        
    """    
    baseline = trmvals[0][1]
    ex = []
    ex.append(baseline)
    extot = 0
    exmain = 0
    if modinfo['type'] == 1:  # Additive ERR TO(1 + T1 + T2 + ...)
        rr = 1
        for i in range(1, len(trmvals)):
            rr += trmvals[i][1]
            termex = baseline * trmvals[i][1]
            extot = extot + termex
            exjnt = 0 * extot
            ex.append(termex)
    if modinfo['type'] == 5:  # Multiplicative ERR TO(1+T1)(1+T2)...
        rr = 1
        for i in range(1, len(trmvals)):
            err = 1 + trmvals[i][1]
            rr *= err
            termex = baseline * trmvals[i][1]
            exmain += termex
            ex.append(termex)
        extot = baseline * ( rr - 1)
        exjnt = extot - exmain           
    if ((modinfo['type'] == 3) or (modinfo['type'] == 2)):  # EAR or Additive RR T0*(T1+T2+...)
        rr = 0
        for i in range(1, len(trmvals)):
            rr += trmvals[i][1]
            if modinfo['type'] == 3:
                ex.append(rr)
                extot += rr
        exjnt = 0 * extot

    if modinfo['type'] == 3: 
        fv = baseline + rr  # EAR
    else:
        fv = baseline * rr  # RR
# returns PY specific total fitted cases, total excess, subterm-specific excess (main effect)
    fv = np.reshape(fv,(-1,1))
    fvinfo = {}
    fvinfo.update({'fv': fv,
        'extot': extot,
        'termfv': ex,
        'exjnt': exjnt})

    return(fvinfo)

# end of calc_fvalues
    
def showfittedcases( modelinfo, obscases, fvinfo, tnames= ""):
    # Need to update for different model types
    nterms = len(fvinfo['termfv'])
    if tnames == "":
        tnames= "Baseline"
        for tno in range(1,nterms):
            tn = "Term"+str(tno)
            tnames = np.append(tnames, tn)
    print '\n', '%36s' % 'Observed cases:', '%7.0f' % sum(obscases)
    print '%36s' % '  Fitted cases:', '%10.2f' % np.sum(fvinfo['fv'])
    for i in range(len(tnames)):
        print '%22s' % tnames[i], "fitted cases: ", '%10.2f' % sum(fvinfo['termfv'][i])
    


    