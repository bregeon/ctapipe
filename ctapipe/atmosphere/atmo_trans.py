#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read the atmospheric transmission and extinction coefficient
for a given Modtran like atmospheric transmission table

@author: Johan Bregeon
"""
import logging
logger = logging.getLogger(__name__)

import numpy
from scipy import interpolate


def readAtmoTrans(file_path='../data/atm_trans_1800_1_10_0_0_1800.dat'):
    ''' Read atmospheric transmission table

    Parameters
    ----------
    file_path : string
        atmospheric profile file name
    '''
    content=open(file_path).readlines()

    opt_depth = {}
    extinction = {}
    for l in content[14:]:
        all = l.split()
        if l[0] == '#':
            obs_alt = float(all[2].strip(','))
            altitudes = numpy.array(all[4:],'f') #a.s.l.
            altitudes = altitudes*1000 # in meters
        else:
            opt_depth[float(all[0])] = numpy.array(all[1:],'f')

    for wl,depth in opt_depth.items():
        alpha=[]
        for i in range(len(depth)-1):
            alpha.append((depth[i+1]-depth[i])/(altitudes[i+1]-altitudes[i]))
        extinction[wl] = numpy.array(alpha,'f')

    return content[:15], obs_alt, altitudes, opt_depth, extinction

def writeAtmoTrans(run, header, scaled_od, fname=None):
    ''' Dump the scaled optical depth table to an ASCII file
    in Kaskade/Corsika input file format

    Parameters
    ----------
    run : int
        a number to identify the correspoding data run
    header : string
        the header of modtran like atmospheric transmission files
    scaled_od : dict
        dictionary of list of optical depth vs height per wavelength
    fname : string
        explicitely give the output file name
    '''
    if fname is None:
        fname = 'atm_trans_Lidar_%s.dat'%run
    content = header
    for wl, od_prof in scaled_od.items():
        txt='        %03d     '%wl
        for od in od_prof:
            if od<10:
                txt += '%10.6f'%od
            else:
                txt += '%10.2f'%od
        content.append(txt+'\n')
    open(fname,'w').writelines(content)
    logger.info('%s written.' % fname)

def getODForAltitudes(npalt, npod, newalt, opt_depth_wl):
    '''
    Get Transmission values for a given array of altitudes

    Parameters
    ----------
    npalt : array
        array of altitudes
    npod : array
        array of optical depths
    newalt : array
        array of interpolated altitudes
    opt_depth_wl : array
        array of optical depth for a given wavelength from the model
    '''
    # x=alt y=nptrans
    tck = interpolate.splrep(npalt, npod, s=0)
    newod = interpolate.splev(newalt, tck, der=0)
    # no extrapolation, use model outside measurements
    deltaod = 0
    for i in range(len(newalt)):
        #print('%s alt %s %s %s'%(i, newalt[i], npalt[0], npalt[-1]))
        if newalt[i] < npalt[0]: #below AltMin - constant
            newod[i] = npod[0]
        elif newalt[i] > npalt[-2]: #above R0 -1 bin - follow model
            if deltaod == 0:
                deltaod = newod[i]-opt_depth_wl[i]
                #print('deltaod %s %s %s'%(deltaod, i, newalt[i]))
            newod[i] = opt_depth_wl[i]+deltaod
    return newod

def getODRatio(newod, opt_depth_wl):
    """
    calculate the ratio of optical depth profile for a given wavelength

    Parameters
    ----------
    newod : array
        array of optical depth for a given wavelength
    opt_depth_wl : array
        array of optical depth for a given wavelength from the model
    """
    od_ratio = newod/opt_depth_wl
    return od_ratio

def getScaledOD(opt_depth, od_ratio, corr=1):
    '''
    Scale all wavelength using 532 nm measurement

    Parameters
    ----------
    opt_depth : dict
        dictionary of optical depth profile per wavelength
    od_ratio : array
        array of ratio of data over model to be applied for correction
    '''
    scaled_od = {}
    for wl,od in opt_depth.items():
        scaled_od[wl] = opt_depth[wl]*od_ratio*corr #pow(wl/532.,1.3)
        for i in range(len(scaled_od[wl])):
            if scaled_od[wl][i] > 10:
                scaled_od[wl][i] = 99999.00
    return scaled_od

def getIndex(myarray,val):
    '''
    finds the index of first value above a given one in a profile

    Parameters
    ----------
    myarray : array
        array of monotonic values, typically an atmospheric profile
    val : float
        threshold value
    '''
    return len(myarray)-sum(myarray>val)

if __name__=="__main__":
    ''' Write up an example of reading
    '''
