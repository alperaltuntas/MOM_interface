#!/usr/bin/env python3

import os
import xarray as xr
import numpy as np
import argparse

descr = """
Utility functions for MOM6 land block elimination preprocessor module
"""

def MOM_define_layout(isz, jsz, ndivs):
    """ This function is a Python implementation of MOM_define_layout subroutine from 
    MOM_domains.F90, and is used here for guiding the land block elimination functions.
    Given a global array size (isz x jsz) and a number of (logical) processors (ndivs)
    , provide a layout of the processors in the two directions where the total number
    of processors is the product of the two layouts and number of points in the 
    partitioned arrays are as close as possible to an aspect ratio of 1. """

    # First try to divide ndivs to match the domain aspect ratio.  If this is not an even
    # divisor of ndivs, reduce idiv until a factor is found.
    idiv = max(np.rint(np.sqrt( (ndivs*isz)/jsz) ), 1)
    while ndivs%idiv != 0:
        idiv -= 1
    jdiv = ndivs // idiv

    return int(idiv), int(jdiv) 

def mpp_compute_extent(isg, ieg, ndivs):
    """ This function is a Python implementation of mpp_compute_extent from FMS 
    mpp_domains_define.inc and is used by builnml for generating mask tables for
    land block elimination.
    Computes extents for a grid decomposition with the given indices and divisions"""

    def even(x):
        assert isinstance(x,int)
        return x%2==0
    def odd(x):
        return not even(x)
    
    ibegin = [None for i in range(ndivs)]
    iend = [None for i in range(ndivs)]
    
    is_ = isg

    symmetrize = ( even(ndivs) and even(ieg-isg+1) ) or \
        ( odd(ndivs) and odd(ieg-isg+1) ) or \
        ( odd(ndivs) and even(ieg-isg+1) and ndivs<(ieg-isg+1)/2 )
        
    imax = ieg
    ndmax = ndivs
    
    for ndiv in range(ndivs):

        # do bottom half of decomposition, going over the midpoint for odd ndivs
        if ndiv<(ndivs-1)//2+1:
            ie = is_ + int(np.ceil( (imax-is_+1)/(ndmax-ndiv) )) - 1
            ndmirror = (ndivs-1) - ndiv # mirror domain
            if ndmirror > ndiv and symmetrize:
                # mirror extents, the max(,) is to eliminate overlaps
                ibegin[ndmirror] = max(isg+ieg-ie, ie+1)
                iend[ndmirror] = max(isg+ieg-is_, ie+1)
                imax = ibegin[ndmirror] - 1
                ndmax -= 1
        else:
            if symmetrize:
                # do top half of decomposition by retrieving saved values
                is_ = ibegin[ndiv]
                ie = iend[ndiv]
            else:
                ie = is_ + int(np.ceil( (imax-is_+1)/(ndmax-ndiv) )) - 1
        
        ibegin[ndiv] = is_
        iend[ndiv] = ie
        
        assert ie>=is_, "domain extents must be positive definite.'"
        assert not(ndiv==ndivs-1 and iend[ndiv]!=ieg), "domain extents do not span space completely"
        
        is_ = ie+1

    assert None not in ibegin, "Error in mpp_compute_extent"
    assert None not in iend, "Error in mpp_compute_extent"
    return ibegin, iend
