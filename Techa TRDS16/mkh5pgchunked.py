# -*- coding: utf-8 -*-
"""
Created on Sun Sfp 10 10:48:27 2017

@author: prest
"""

import os
import sys
import numpy as np   # numerical methods
import h5py  # Hdf5 library
import pandas as pd

sys.path.append('f:/TRDS16MC_HDF5_progs')

import hdf5utils as trch5


dsrlzindir = "f:/trds16mc1709/"
dsrlzh5dir = "f:/hdf5/"

# from hdf5utils import addtoh5pg


os.chdir(dsrlzh5dir)



# batchstats = trch5.batchinfo(icsvfn, iext)



# TRC HP5 person group file
batchpfx = "2017Cases"
batchids = [ '1-1000', '1k-2k'  ,   '2k-3k', '3k-4k'  , '4k-5k'  , 
              '5k-6k', '6k-7k'  ,   '7k-8k', '8k-9k'  , '9K-10K' ,  
            '10k-11k', '11K-12K', '12k-13k', '13k-14k', '14k-15k',
            '15k-16k', '16k-17k', '17k-18k', '18k-19k', '19k-20k',
            '20k-21k', '21k-22k', '22k-23k', '23k-24k', '24k-25k', 
            '25k-26k', '26k-27k', '27k-28k', '28k-29k', 
            'makeups']


cohort = "trc"
h5filenm = dsrlzh5dir + cohort + '16mc.h5'
batchlist = ['1-1000', '1k-2k'  ,   '2k-3k', '3k-4k'  , '4k-5k'  , 
              '5k-6k', '6k-7k'  ,   '7k-8k', '8k-9k'  , '9K-10K' ,  
            '10k-11k', '11K-12K', '12k-13k', '13k-14k', '14k-15k',
            '15k-16k', '16k-17k', '17k-18k', '18k-19k', '19k-20k',
            '20k-21k', '21k-22k', '22k-23k', '23k-24k', '24k-25k', 
            '25k-26k', '26k-27k', '27k-28k', '28k-29k', '29k-30k', 
            'makeups']
for batchid in batchlist:
    batchnm  = dsrlzindir + batchpfx + batchid 
    print batchid
    trch5.h5update(h5filenm,cohort,batchid,batchnm) 

'''
cohort = "eurt"
h5filenm = dsrlzh5dir + '/x' + cohort + '16mc.h5'
# batchlist = ['26k-27k', '27k-28k', '28k-29k']
for batchid in batchlist:
    batchnm  = dsrlzindir + batchpfx + batchid 
    print batchid
    trch5update(h5filenm,cohort,batchid,batchnm)
'''

